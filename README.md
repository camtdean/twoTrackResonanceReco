# twoTrackResonanceReco

Reconstruction of simple two-body resonance decays directly from already-fitted sPHENIX tracks without using Kalman filters (though I don't know if the ACTS propagator uses one, possibly?). It pairs opposite-sign tracks, finds their secondary vertex with the ACTS track propagator, applies kinematic selections, and writes a ROOT nTuple of candidates.

## Requirements

- A DST (real data or simulation) with `SvtxTrackMap` and `SvtxVertexMap` on
  the node tree.
- `ActsGeometry` **must** be on the node tree. This is a hard requirement —
  the module aborts the event if it can't find it, because every DCA and
  vertex calculation in this module goes through the ACTS propagator. This
  means your macro needs to load the tracking geometry
  (`G4_ActsGeom.C` / `TrackingInit()`) even if you're only reading an
  already-reconstructed DST and never touch raw hits.
- `TRKR_CLUSTER` (or `TRKR_CLUSTER_SEED`) and `TPCGEOMCONTAINER` are optional
  — they're only needed if you turn on dE/dx-based PID (`usePID(true)`). If
  they're missing, the module just skips dE/dx and keeps running.

## Building

This follows the usual autotools pattern for an sPHENIX analysis package:

```
./autogen.sh
mkdir build && cd build
../configure --prefix=$MYINSTALL
make -j $(nproc) install
```

## How it works

Each event, `process_event()` does the following, in order:

1. **Node retrieval.** `SvtxTrackMap`, `SvtxVertexMap`, and `ActsGeometry`
   are fetched; any one missing aborts the event. The cluster container and
   TPC geometry container are fetched too, but only to enable dE/dx — their
   absence just turns PID off for the run.

2. **Track selection.** Every track in `SvtxTrackMap` is required to:
   - have at least one MVTX or INTT cluster state — this is because the original intent is to reconstruct D0
   - pass `chi2/ndf < ` the value set by `setMaxDaughterChi2perNDF()` (default is 100);
   - pass `pT > ` the value set by `setMinDaughterPT()` (default is 0).

3. **Pairing.** Every combination of two `goodTracks` is considered. A pair
   is rejected immediately if the two tracks have the same charge sign, or
   if they weren't reconstructed on the same beam crossing. For pairs that
   survive, the module looks up the primary vertex on that same crossing
   (picking the one with the best `chi2/ndof` if there's more than one).

4. **Secondary vertex finding (`buildSV`).** Rather than a covariance-weighted vertex fit, it's an
   iterative fixed-point search built entirely on ACTS propagation:
   - seed a target point at the midpoint of the two tracks' raw
     `(x, y, z)` positions;
   - propagate both tracks to that target point with the ACTS propagator
     (`makeTrackParams` + `makeVertexSurface` + `propagateTrackFast`);
   - move the target to the midpoint of the two propagated positions;
   - repeat (up to 10 iterations, or until the target moves less than
     `1e-4` cm between iterations).

   The result is a secondary vertex position, each track's momentum at that
   point, and the track-to-track DCA (the separation between the two
   propagated positions at convergence).

5. **Track-to-track DCA cut** — `setDaughterDCACut()`.

6. **Flight distance cut** — the SV-to-PV separation, via
   `setFlightDistanceCut()`.

7. **DIRA cut** — cosine of the angle between the reconstructed mother
   momentum and the PV→SV flight direction, via `setDIRACut()`.

8. **Mother PV impact parameter cut** — how far the mother's flight line
   misses the primary vertex, via `setMotherIPCut()`.

9. **Daughter PV DCA cut.** Each daughter is separately propagated back to
   the primary vertex (`trackToVertexDCA`), and the *smaller* of the two
   daughter DCAs must clear `setDaughterIPCut()`. 

10. **Mass hypothesis assignment.** If the two daughter species are
    identical (e.g. K-short → π+π-), the charge sign alone fixes which
    track is "daughter 1" vs "daughter 2." If they're different species
    (e.g. D0 → K-π+), both charge assignments are tried as independent
    hypotheses, each gets its own invariant mass, and each is checked
    against `setMotherMassRange()` on its own.

11. **Optional dE/dx PID** — `usePID(true)`. Measured dE/dx
    (`TrackAnalysisUtils::calc_dedx`, using TPC clusters) is compared
    against expected dE/dx bands loaded from the `TPC_DEDX_FITPARAM`
    calibration payload, separately for pion/kaon/proton and each charge
    sign. A candidate fails PID if either daughter's measured dE/dx falls
    outside `setdEdxBandWidth()` (a fractional width, default 20%) of the
    expected band for its assumed species.

12. **Tree fill.** Every hypothesis that survives all of the above gets its
    own row in the output tree. For non-identical daughter species where
    *both* charge hypotheses pass, `both_charge_states_passed` is set to
    `true` on both rows so you can identify and handle that ambiguity
    downstream.

## Quick start

```cpp
twoTrackResonanceReco *myKshortReco = new twoTrackResonanceReco("KshortReco");
myKshortReco->setMotherMassRange(0.4, 0.6);   // GeV, around the K-short mass
myKshortReco->setDaughterDCACut(0.1);          // cm
myKshortReco->setFlightDistanceCut(0.8);       // cm
myKshortReco->setOutputFileName(outputFileName);
se->registerSubsystem(myKshortReco);
```

Daughter species default to `211, -211` (π+π-), so the snippet above is
already a working K-short reconstruction. For a different decay, set the
daughter PDG IDs and adjust the cuts to match, for example:

```cpp
twoTrackResonanceReco *myD0Reco = new twoTrackResonanceReco("D0Reco");
myD0Reco->setDaughterPDGIDs(321, -211);        // K- , pi+  (and the charge conjugate)
myD0Reco->setMotherMassRange(1.7, 2.0);        // GeV, around the D0 mass
myD0Reco->setDaughterDCACut(0.05);
myD0Reco->setFlightDistanceCut(0.02);
myD0Reco->usePID(true);
se->registerSubsystem(myD0Reco);
```

A full, runnable example — including the mandatory tracking-geometry setup
— lives in `macro/Fun4All_twoTrackReco.C`.

## Configuration reference

| Setter | Meaning | Default |
|---|---|---|
| `setDaughterPDGIDs(pdgID1, pdgID2)` | PDG IDs of the two daughter species (sign matters — fixes which charge is "daughter 1") | `211, -211` |
| `setMotherMassRange(min, max)` | Invariant mass window a candidate must fall in [GeV] | `0, 2` |
| `setMaxDaughterChi2perNDF(cut)` | Max track fit `chi2/ndf` for a daughter to be considered | `100` |
| `setMinDaughterPT(cut)` | Min daughter track `pT` [GeV] | `0` (off) |
| `setDaughterDCACut(cut)` | Max track-to-track DCA at the secondary vertex [cm] | `999` (off) |
| `setFlightDistanceCut(cut)` | Min PV-to-SV separation [cm] | `-999` (off) |
| `setMotherIPCut(cut)` | Max mother impact parameter to the PV [cm] | `999` (off) |
| `setDaughterIPCut(cut)` | Min of the two daughters' DCA to the PV [cm] | `-0.1` (off) |
| `setDIRACut(cut)` | Min cosine of the angle between mother momentum and flight direction | `-1.1` (off) |
| `usePID(bool)` | Turn on dE/dx-based PID | `false` |
| `setdEdxBandWidth(width)` | Fractional half-width of the accepted dE/dx band | `0.2` |
| `setTrackMapName(name)` | Node name for the track map | `"SvtxTrackMap"` |
| `setVertexMapName(name)` | Node name for the vertex map | `"SvtxVertexMap"` |
| `setOutputFileName(name)` | Output ROOT file name | `"twoTrackResonanceReco.root"` |

Most of the geometric cuts default to effectively "off" (a value the real
observable can never fail), so you'll want to set the ones relevant to your
decay explicitly rather than relying on the defaults for anything but the
mass window and quality/pT cuts.

## Output tree (`DecayTree`)

One row per candidate that passes every cut.

| Branch | Meaning |
|---|---|
| `crossing` | Beam crossing shared by the two daughter tracks |
| `PV_x, PV_y, PV_z` | Primary vertex position used for this candidate [cm] |
| `SV_x, SV_y, SV_z` | Reconstructed secondary (decay) vertex [cm] |
| `mother_mass` | Reconstructed invariant mass [GeV] |
| `mother_pT, mother_eta, mother_phi` | Reconstructed mother kinematics |
| `mother_DIRA` | Cosine of the pointing angle |
| `mother_flight_distance` | PV-to-SV distance [cm] |
| `mother_PV_DCA` | Mother impact parameter to the PV [cm] |
| `daughter{1,2}_mass` | Mass assumed for that daughter hypothesis [GeV] |
| `daughter{1,2}_charge` | Daughter track charge |
| `daughter{1,2}_pT` | Daughter transverse momentum from the SvtxTrack object |
| `daughter{1,2}_eta` | Daughter pseudorapidity from the SvtxTrack object |
| `daughter{1,2}_phi` | Daughter track phi from the SvtxTrack object |
| `daughter{1,2}_PV_DCA` | That daughter's DCA to the primary vertex [cm] |
| `daughter{1,2}_dEdx` | Measured TPC dE/dx (`-1` if unavailable) |
| `daughter{1,2}_chi2_per_ndf` | That daughter's track fit `chi2/ndf` |
| `track_to_track_DCA` | Daughter-to-daughter DCA at the secondary vertex [cm] |
| `both_charge_states_passed` | `true` if both charge hypotheses passed for a non-identical-species pair |

## Relationship to KFParticle_sPHENIX

If you've used `KFParticle_sPHENIX`, the inputs and general shape of the
selection (mass window, DCA cuts, flight distance, DIRA, PID) will look
familiar — that's intentional. The important difference is underneath:
`KFParticle_sPHENIX` runs an actual Kalman-filtered vertex fit, which
refines each daughter's momentum and covariance at the vertex and produces
genuine chi2-based significances (decay length significance, IP
significance, vertex chi2/ndof). `twoTrackResonanceReco` does not do a
covariance-weighted fit at all — `buildSV()` is a purely geometric,
iterative point-of-closest-approach search using ACTS propagation, with no
covariance propagation or momentum refinement. Every DCA, IP, and flight
distance in this module's output is a raw geometric distance in cm, not a
statistical significance, and there's no vertex chi2/ndof to cut on because
no vertex fit is ever performed.

In practice this means: cut values that were tuned as "N standard
deviations" in a KFParticle analysis are not directly portable here — you're
choosing plain distance cuts instead, and you should expect somewhat looser
background rejection at a given signal efficiency, particularly for cuts
where KFParticle would use a significance. 

## Simulation macro (`simulation_macro/Fun4All_D0_sim.C`)

```cpp
Fun4All_D0_sim(const int nEvents = 10, const std::string &outdir = "./",
               const int processID = 0, bool doPolytracking = false)
```

- `processID` numbers a single simulation job. It's zero-padded to 5 digits
  and baked into every output filename, so file names stay aligned once
  you're running more than 9 jobs in parallel (`_00003.root` rather than
  `_3.root`).
- `doPolytracking` switches the TPC track-reconstruction path:
  - `false` (default): the standard CA-seeding path — `TpcClusterizer` with
    `SetDeadChannelMapName("TPC_DEADCHANNELMAP")` (picks up the real TPC
    dead-channel map for the configured run number, if one exists) followed
    by `Tracking_Reco()`.
  - `true`: the full polytracking/polyseeding chain — silicon seeding, TPC
    module/assembled tracks, `TpcCrossingFinder`, `Tpc_PolyClusterizer`,
    poly track reconstruction and vertexing, then track matching. Both
    `TpcCrossingFinder` and `Tpc_PolyClusterizer` are configured with
    `setIsNonDistortedMC(true)`, which turns off the ExB drift effect the
    polyline drift model would otherwise apply even for undistorted MC.
    **This flag doesn't exist in mainline `coresoftware` yet** — it's from
    a draft PR (sPHENIX-Collaboration/coresoftware#4455), so you need
    `offline/packages/PHGarfield` and `offline/packages/tpctrackreco` built
    from that branch before this macro will run.
- Output from every reconstruction channel registered in the macro
  (currently just `DzeroReco`) is written under a single `output/`
  directory, with a per-channel subdirectory inside it (e.g.
  `output/Dzero_reco_caseeding_/`, `output/Dzero_reco_polyseeding_/`).
  Keeping everything under one `output/` parent means the batch scripts
  below can copy results back from scratch space without needing to know
  the names of individual reconstruction channels.

## Real-data macro (`macro/Fun4All_twoTrackReco.C`)

```cpp
Fun4All_twoTrackReco(const int nEvents = 1000,
                      const std::string &inputList = "jobLists/run79516_00.txt", const int nSkip = 0)
```

Reads a list of DST files (one per line, `inputList`), derives the run
number and segment from the first line to build a zero-padded output
filename, and registers a K-short reconstruction (`KshortReco`) as a
starting example — see Quick start above for how to add more channels.
Like the simulation macro, output goes under a single `output/` directory
(`output/Kshort_reco/`), for the same scratch-copy-back reason. 

## Batch submission with Condor

Two matched pairs of scripts handle batch submission, one for real data and
one for simulation:

| | Real data | Simulation |
|---|---|---|
| Shell wrapper | `macro/runData.sh` | `simulation_macro/runSims.sh` |
| Condor submit file | `macro/submitData.job` | `simulation_macro/submitSim.job` |

Both shell wrappers follow the same pattern: source the sPHENIX
environment, point `LD_LIBRARY_PATH`/`ROOT_INCLUDE_PATH` at your local
install (`$MYINSTALL`), then invoke the matching macro with `root.exe -q -b`.

**`useScratch`.** Each wrapper has a `useScratch` flag (`false` by
default). When `true`, the script rsyncs its own directory into the job's
private `$_CONDOR_SCRATCH_DIR` before running, and rsyncs the `output/`
directory back out to the original shared directory afterward. This keeps
large batch submissions from having every job read and write the same
shared directory concurrently. With `useScratch=false`, everything just
runs directly in `initialDir` (the shared directory named in the submit
file) and nothing needs to be copied back.

**Memory retries.** Both submit files use job retry
patterns to avoid overconsumption of memory: `request_memory` sets the initial allocation, and
`retry_request_memory_increase`/`retry_request_memory_max` let Condor
automatically retry a job with more memory (up to the max) instead of
evicting it outright. 

**Wiring up the job list.**
- `submitData.job` queues one job per line of `macro/inputDSTlists.txt`,
  where each line is itself a path to a file listing up to 100 DST files
  (`macro/jobLists/run<runnumber>_NN.txt`). Splitting a full run's file
  list into chunks like this can be done with:
  `gsplit -l 100 -d --additional-suffix=.txt run<runnumber>.list run<runnumber>_`
  (drop the leading `g` on Linux).
- `submitSim.job` queues jobs by process number (`Queue 1000` submits
  processes `0`–`999`), and each process number becomes the `processID`
  argument passed through to `Fun4All_D0_sim.C`.
