# ---------------------------------------------------------------------------------------#
# takes three3 files:
#       1) dllee_vertex variables
#       2) tracker_reco tree
#       3) 
# ---------------------------------------------------------------------------------------#
import ROOT
from sys import argv
from ROOT import gDirectory, larlite, TFile, TCanvas, gROOT, TH1D, TGraph, TF1, TGraphErrors, TLine, TFunction, Double
from ROOT import TH3F, TH3D
import numpy as np
import spatial_cal

MinimumTrackLength = 5
binsize = 10
view = 2
def InVolume(x,y,z,edge_x=10,edge_y=10,edge_z=10):

    # ---------------------------------------------------------------------------------------#
    # Function uses detector active volume coordinates in conjuction with an input parameter
    # for distance from edge of active volume to define a fiducial volume and return 0/1 if
    # the given coordinates are contained in the so defined fiducial volume, default 10 cm
    # ---------------------------------------------------------------------------------------#       

    #--Defines active volume coordinate boundaries--#                                       
    xmax =  256.25
    xmin =  0
    ymax =  116.5
    ymin = -116.5
    zmax =  1036.8
    zmin =  0
    #-----------------------------------------------#          

    if x < (xmax - edge_x) and x > (xmin + edge_x) and y < (ymax - edge_y) and y > (ymin + edge_y) and z < (zmax - edge_z) and z > (zmin + edge_z):
        return 1
    else:
        return 0

# --- load tree with cut information and stores info in a [run, subrun, event] vector for reference --- #
root_file = argv[1]
myfile = TFile(root_file)

mychain = gDirectory.Get("NuMuVertexVariables")
entries = mychain.GetEntriesFast()
runs, good3drun, cosmicll, runsub = [], [], [], []

for entry in range( entries ):
    ientry = mychain.LoadTree( entry ) 
    if ientry < 0:
        break

    nb = mychain.GetEntry( entry )
    if nb <= 0:
        continue

    run  = mychain.run
    subrun = mychain.subrun
    event = mychain.event
    vtxid = mychain.vtxid
    goodrun = mychain.Good3DReco
    cosll = mychain.CosmicLL
    run_vector = [run, subrun, event, vtxid]
    good3drun.append(goodrun)
    cosmicll.append(cosll)
    runs.append(run_vector)

runs = np.asarray(runs)

# --- initializing some arrays --- #
thrumu_rr, thrumu_x, thrumu_y, thrumu_z, thrumu_dQdx, thrumu_length = [], [], [], [], [], []
protrr, protx, proty, protz, protdQdxY = [], [], [], [], []

# ---------------------------------------------------------------------------------------#
# loads track_reco tree and selects the crossing comsmic tracks (thrumu) for the spatial
# calibration and proton tracks for the recombination parameters.
# thrumu -- selects based on:
#           1) length of combined tracks about size of the detector
#           2) length of individual tracks greater than 30 cm (to account for Michel's)
# protons -- selects based on:
#           1) 2 tracks from vertex, chooses higher, average dQdx AND shorter length
#           2) good3dreco
#           3) cosmicll > 4
# ---------------------------------------------------------------------------------------#

track_infile = argv[2]
track_manager = larlite.storage_manager()
track_manager.set_io_mode(track_manager.kREAD)
track_manager.add_in_filename(track_infile)
track_manager.set_in_rootdir("")
track_manager.open()

while track_manager.next_event():
    track = track_manager.get_data(larlite.data.kTrack,"trackReco")
    vertex = track_manager.get_data(larlite.data.kVertex,"trackReco")
    if track.size() != 0L:
        vert = []
        xprev = track[0].Vertex().x()
        vertex_number = 0
        for idx in range(len(track)):
            xver = track[idx].Vertex().x()
            if xver == xprev:
                vert.append(vertex_number)
            if xver != xprev:
                vertex_number += 1
                vert.append(vertex_number)
                xprev = xver

        thrumu_track, protontrack, muontrack = [], [], []
        for idx in range(len(vertex)):
            cosmic_inttrack = []
            for i in range(len(track)):
                if vert[i] == idx:
                    length = track[i].Length()
                    if length >= 30:
                        cosmic_inttrack.append(i)
            if len(cosmic_inttrack) != 2:
                break
            Length = []
            Y_pos = []
            for i in cosmic_inttrack:
                length = track[i].Length()
                t = int(track[i].NumberTrajectoryPoints()) - 1
                y_pos = track[i].LocationAtPoint(t).Y()
                Length.append(length)
                Y_pos.append(y_pos)
            L = Length[0] + Length[1]
            if L > 225 and L < 280 and Length[0] > 35 and Length[1] > 35:
                thrumu_track.append(cosmic_inttrack[0])
                thrumu_track.append(cosmic_inttrack[1])
                thrumu_length.append(L)
        for idx in range(len(track)):
            rr, dQdx = [], []
            trk = track[idx]
            for i in range(int(trk.NumberTrajectoryPoints())):
                rr.append(trk.Length(i))
                dQdx.append(2.1*trk.DQdxAtPoint(i,view))
            if idx in thrumu_track:
                for f in range(len(rr)):
                    thrumu_x.append(trk.LocationAtPoint(f).X())
                    thrumu_y.append(trk.LocationAtPoint(f).Y())
                    thrumu_z.append(trk.LocationAtPoint(f).Z())
                thrumu_rr.extend(rr)
                thrumu_dQdx.extend(dQdx)

    if track.size() != 0L:
        vert = []
        xprev = track[0].Vertex().x()
        vertnum = 0
        for idx in range(len(track)):
            xver = track[idx].Vertex().x()
            if xver == xprev:
                vert.append(vertnum)
            if xver != xprev:
                vertnum = vertnum + 1
                vert.append(vertnum)
            xprev = xver

        # ------ Assign Each Track as Either a Proton or Muon ------ #

        for idx in range(len(vertex)):
            proton_inttrack = []
            for i in range(len(track)):
                if vert[i] == idx:
                    vtx_x = track[i].Vertex().x()
                    vtx_y = track[i].Vertex().y()
                    vtx_z = track[i].Vertex().z()
                    vtxInFiducial = InVolume(vtx_x, vtx_y, vtx_z)
                    length = track[i].Length()
                    if length >= MinimumTrackLength and vtxInFiducial ==True:
                        proton_inttrack.append(i)
            if len(proton_inttrack) != 2: # Only keep vertexes with two tracks
                break
            # --- Find Average dQdx of Each Track --- #
            AdQdx = []
            Tot_L2 = []
            for i in proton_inttrack:
                trk = track[i]
                points = int(trk.NumberTrajectoryPoints())
                deadpoints, avgdQdx = 0, 0
                Tot_L = trk.Length()
                Tot_L2.append(Tot_L)
                for f in range(points):
                    for g in range(3):
                        dQdxpoint = trk.DQdxAtPoint(f,g)
                        avgdQdx += dQdxpoint
                        if dQdxpoint == 0:
                            deadpoints += 1
                AdQdx.append(avgdQdx/(((3*points) - deadpoints)))
            # --- Assign Higher Average dQdx Tracks to Proton, Lower to Muon --- #
            if AdQdx[0] > AdQdx[1]:
                if Tot_L2[0] < Tot_L2[1]:
                    protontrack.append(proton_inttrack[0])
            else:
                if Tot_L2[1] < Tot_L2[0]:
                    protontrack.append(proton_inttrack[1])

        # ------ For Each Proton/Muon Track, Append dQdx, Location, and Length ------ #
        for idx in range(len(track)):
            pos, dQdxU, dQdxV, dQdxY = [], [], [], []
            trk = track[idx]
            num = int(track[idx].NumberTrajectoryPoints())
            # --- Check to See How Many Deadwires are in a Track, Only Keep Tracks Less Than Some Number of Deadwires --- #
            eventid = int(track.event_id())
            run = int(track.run())
            subrun = int(track.subrun())
            vertexid = int(vertex[vert[idx]].ID())
            run_vector = [run, subrun, eventid, vertexid]
            run_vector = np.asarray(run_vector)
            runsub.append([run, subrun])
            if run_vector in runs:
                runloc = np.logical_and(run_vector[0] == runs[:,0], run_vector[1] == runs[:,1])
                subrunloc = np.logical_and(run_vector[2] == runs[:,2], run_vector[3] == runs[:,3])
                ind = np.logical_and(runloc == True, subrunloc == True)
                ibin = np.where(subrunloc == True)[0][0]
                if idx in protontrack:
                    if good3drun[ibin] == 1 and cosmicll[ibin] > 4: 
                        for v in range(num):
                            f_len = trk.Length()
                            c_len = trk.Length(v)
                            cur_len = f_len - c_len
                            if cur_len > 2:
                                protx.append(trk.LocationAtPoint(v).X())
                                proty.append(trk.LocationAtPoint(v).Y())
                                protz.append(trk.LocationAtPoint(v).Z())
                                protrr.append(trk.Length(v))
                                protdQdxY.append(trk.DQdxAtPoint(v, view))


track_manager.close()

# ---------------------------------------------------------------------------------------#
# spatial calibration - removes the deadwire hits from the thrumu tracks and then finds
# calibration factors for YZ plane and X axis.  then applies calibration factors to 
# proton tracks
# ((uses spatial_cal function))
# ---------------------------------------------------------------------------------------#

protx, proty, protz, protdQdxY, protrr = np.asarray(protx), np.asarray(proty), np.asarray(protz), np.asarray(protdQdxY), np.asarray(protrr)
thrumu_rr, thrumu_x, thrumu_y, thrumu_z, thrumu_dQdx = np.asarray(thrumu_rr), np.asarray(thrumu_x), np.asarray(thrumu_y), np.asarray(thrumu_z), np.asarray(thrumu_dQdx)

ind = np.where(thrumu_dQdx != 0)
thrumu_rr, thrumu_x, thrumu_y, thrumu_z, thrumu_dQdx = thrumu_rr[ind], thrumu_x[ind], thrumu_y[ind], thrumu_z[ind], thrumu_dQdx[ind]
ind = np.where(protdQdxY != 0)
protrr, protx, proty, protz, protdQdxY = protrr[ind], protx[ind], proty[ind], protz[ind], protdQdxY[ind]


grid_zy, bins, binloc, zmesh, ymesh = spatial_cal.yz_data(thrumu_z, thrumu_y, thrumu_dQdx, binsize=binsize)

globalval = 100
localval_zy = np.divide(globalval, grid_zy)


thrumu_dQdx2 = np.copy(thrumu_dQdx)
nrow, ncol = grid_zy.shape
for row in range(nrow):
    for col in range(ncol):
        zc = zmesh[row, col]
        yc = ymesh[row, col]

        posz = np.abs(thrumu_z - zc)
        posy = np.abs(thrumu_y - yc)
        ibin = np.logical_and(posz < binsize/2., posy < binsize/2.)
        ind = np.where(ibin == True)[0]
        calfac = localval_zy[row][col]
        if np.isfinite(calfac) == True:
            thrumu_dQdx2[ind] = thrumu_dQdx[ind]*calfac

grid_x, X = spatial_cal.x_data(thrumu_x, thrumu_dQdx2, binsize=binsize)


localval_x = np.divide(globalval, grid_x)


protdQdxY2 = np.copy(protdQdxY)
nrow, ncol = grid_zy.shape
for row in range(nrow):
    for col in range(ncol):
        zc = zmesh[row, col]
        yc = ymesh[row, col]

        posz = np.abs(protz - zc)
        posy = np.abs(proty - yc)
        ibin = np.logical_and(posz < binsize/2., posy < binsize/2.)
        ind = np.where(ibin == True)[0]
        localcal = localval_zy[row][col]
        if np.isfinite(localcal) == True:
            protdQdxY2[ind] = protdQdxY[ind]*localcal

protdQdxY3 = np.copy(protdQdxY2)
for row in range(len(X)):
    xc = X[row]
    posx = np.abs(protx - xc)
    ibin = np.logical_and(posx < binsize/2., posx < binsize/2.)
    ind = np.where(ibin == True)[0]
    localcal = localval_x[row]
    if np.isfinite(localcal) == True:
        protdQdxY3[ind] = protdQdxY2[ind]*localcal


# ---------------------------------------------------------------------------------------#
# recombination parameters - takes corrected proton tracks and finds the expected dEdx 
# value at each point based on the residual range.  
# fits modified box model to dqdx v dedx curve to find recombination parameters
# ((assumes $DLLEE_UNIFIED_BASEDIR is set))
# ---------------------------------------------------------------------------------------#


ind = np.where(protrr != 0)
protrr, protdQdxY3 = protrr[ind], protdQdxY3[ind]

dEdx = TFile("$DLLEE_UNIFIED_BASEDIR/larcv/app/Reco3D/Proton_Muon_Range_dEdx_LAr_TSplines.root")
treeprot = gROOT.FindObject("sProtonRange2dEdx")

dEdx = []
for i in range(len(protrr)):
    rr = protrr[i]
    dedx = treeprot.Eval(rr)
    dEdx.append(dedx)
dEdx = np.asarray(dEdx)


hist1d = TH1D( "slice", "slice", 20, 0, 300)
Emin, Emax = dEdx.min(), dEdx.max()
E = np.linspace(3, 14, 30)

mu = []
sigma = []
muerror = []
sigmaerror = []
ebin = []
for i in range(len(E) - 5):
    hist1d.Reset()
    preind = np.logical_and(dEdx > E[i], dEdx <= E[i+1])
    bind = (E[i+1] + E[i])/2
    ebin.append(bind)
    ind = np.where(preind == True)
    for idx in protdQdxY3[ind]:
        hist1d.Fill(idx)
    mean = hist1d.GetMean()
    RMS = hist1d.GetRMS()
    minrange = mean - RMS - 15
    maxrange = mean + RMS
    g = TF1('g', 'gaus', 0, 500)
    hist1d.Fit('g', 'R', '', minrange, maxrange)
    mu.append(g.GetParameter(1))
    sigma.append(g.GetParameter(2))
    muerror.append(g.GetParError(1))
    sigmaerror.append(g.GetParError(2))


mu, sigma, ebin = np.asarray(mu), np.asarray(sigma), np.asarray(ebin)
ey = (ebin[1] - ebin[0])/2
n = len(ebin)
errory = np.full((n,1), ey)
g2 = TGraphErrors(ebin.size, ebin.astype(np.double), mu.astype(np.double), errory.astype(np.double), sigma.astype(np.double))


my_func = TF1('my_func', '(1./[2])*(log((([1]/(1.38*0.273))*x) + [0]))/(([1]/(1.38*0.273))*0.0000236)', 0, 20)
my_func.SetParameters(1.,1.,1000.)
my_func.SetParLimits(0, 0.5, 1.0)
my_func.SetParLimits(1, 0.0, 0.5)
my_func.SetParLimits(2, 500, 3000)
my_func.SetParNames("alpha","beta","C")
fit = g2.Fit('my_func')


p0 = my_func.GetParameter(0)
p1 = my_func.GetParameter(1)
p2 = my_func.GetParameter(2)


# ---------------------------------------------------------------------------------------#
# fills a 3d histogram based off the spatial calibrations and then saves the 3dhist and 
# 3 paramater values to calibration file (any empy bins are assigned a value of 1, saved as Cal_run_run_subrun_subrun
# to do: expand to 3 planes
# ---------------------------------------------------------------------------------------#

zbin, xbin, ybin = len(grid_zy[0]), len(grid_x), len(grid_zy)
xmin, xmax = X[0], X[xbin-1]
zmin, zmax = zmesh[0][0], zmesh[0][zbin-1]
ymin, ymax = ymesh[0][0], ymesh[ybin-1][0]


hist3d = TH3D( "Calibration","Calibration",xbin,xmin,xmax,ybin,ymin,ymax,zbin,zmin,zmax)

for row in range(nrow):
    for col in range(ncol):
        for xrow in range(len(X)):
            zc = zmesh[row, col]
            yc = ymesh[row, col]
            xc = X[xrow]
            cal = localval_zy[row, col]
            cal2 = localval_x[xrow]
            Cal = cal*cal2
            if np.isfinite(Cal) == False:
                Cal = 1
            bin_num = hist3d.FindBin(xc, yc, zc)
            cur_val = hist3d.GetBinContent(bin_num)
            if cur_val == 0:
                hist3d.SetBinContent(bin_num, Cal)
            elif cur_val !=  0: 
                print(bin_num)
                print(Cal)
for x in range(len(hist3d)):
    value = hist3d.GetBinContent(x)
    if value == 0.0:
        hist3d.SetBinContent(x, 1)

file_name = "Cal_" + str(runsub[0][0])+'_' + str(runsub[-1][0]) + '_' + str(runsub[0][1]) + '_' + str(runsub[-1][1]) + ".root"
file_name2 = "Cal_" + str(runsub[0][0]) +'-' + str(runsub[-1][0]) + '.root'
dir_name  = "Cal_" + str(runsub[0][0])+'_' + str(runsub[-1][0]) + '_' + str(runsub[0][0]) + '_' + str(runsub[-1][1])
hist3d_base = TH3D( "Calibration","Calibration",xbin,xmin,xmax,ybin,ymin,ymax,zbin,zmin,zmax)


Calibration = TFile(file_name,"UPDATE")
Calibration.mkdir(dir_name)
Calibration.cd(dir_name)
hist3d.Write("SpatialCalibration")

alpha = ROOT.RooDouble(p0)
beta = ROOT.RooDouble(p1)
C = ROOT.RooDouble(p2)
alpha.Write("alpha")
beta.Write("beta")
C.Write("C")

Calibration.Write()





