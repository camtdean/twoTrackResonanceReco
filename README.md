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
make -j4 install
```

## How it works

Each event, `process_event()` does the following, in order:

1. **Node retrieval.** `SvtxTrackMap`, `SvtxVertexMap`, and `ActsGeometry`
   are fetched; any one missing aborts the event. The cluster container and
   TPC geometry container are fetched too, but only to enable dE/dx — their
   absence just turns PID off for the run.

2. **Track selection.** Every track in `SvtxTrackMap` is required to:
   - have at least one MVTX or INTT cluster state — this is because the original intent is as a D0 searcher
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
| `daughter{1,2}_pT, _eta` | Daughter kinematics at the secondary vertex |
| `daughter{1,2}_phi_beamline` | Daughter track phi at its DCA to the beamline |
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