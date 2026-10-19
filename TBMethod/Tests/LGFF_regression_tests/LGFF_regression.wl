(* LGFF complex-energy regression suite.
   Run from a fresh Mathematica kernel:
   Get["/path/to/LGFF_regression.wl"]
   Place this file beside the revised LGFF.wl.
   Expected: all errors below 1.*^-8 (ordinary, well-conditioned tests).
*)
Get[FileNameJoin[{DirectoryName[$InputFileName], "LGFF.wl"}]];

Block[{e, dims, nBlocks, nSites, boundaries, rand, herm,
  hdiag, sdiag, hlower, slower, sigmaLast, zeroBlock, hcsr,
  scsr, sigmaCSR, direct, full, diagonal, lastColumn, lastBlock,
  expectedDiagonal, expectedColumn, ref, gdirect, leadH0,
  leadS0, leadH1, leadS1, couplingH, couplingS, leadEnergy,
  leadG, leadSigma, sigmaReference, finiteLead, orthogonalDirect,
  orthogonalDiagonal, orthogonalColumn, orthogonalLead,
  leadNumCells = 60, report = <||>, record, testError, k, j},

  SeedRandom[20261019];
  rand[r_, c_] := RandomReal[{-1, 1}, {r, c}] +
    I RandomReal[{-1, 1}, {r, c}];
  herm[r_] := With[{a = rand[r, r]},
    (a + ConjugateTranspose[a])/2];
  record[name_String, measured_, tolerance_:1.*^-8] := (
    AssociateTo[report,
      name -> <|"Error" -> N[measured],
        "Pass" -> TrueQ[measured < tolerance]|>];
  );

  (* 1. Three-block CSR in a nonorthogonal basis, at complex energy. *)
  dims = {2, 3, 2};
  nBlocks = Length[dims];
  nSites = Total[dims];
  boundaries = FoldList[Plus, 0, dims];
  e = 0.17 + 0.035 I;
  hdiag = 0.18 herm[#] & /@ dims;
  sdiag = (IdentityMatrix[#] + 0.012 herm[#]) & /@ dims;
  hlower = Table[0.045 rand[dims[[j + 1]], dims[[j]]],
    {j, nBlocks - 1}];
  slower = Table[0.012 rand[dims[[j + 1]], dims[[j]]],
    {j, nBlocks - 1}];
  sigmaLast = -0.08 I IdentityMatrix[Last[dims]];
  zeroBlock[i_, j_] := ConstantArray[0., {dims[[i]], dims[[j]]}];

  hcsr = ArrayFlatten[Table[
    Which[i == j, hdiag[[i]], i == j + 1, hlower[[j]],
      j == i + 1, ConjugateTranspose[hlower[[i]]],
      True, zeroBlock[i, j]],
    {i, nBlocks}, {j, nBlocks}]];
  scsr = ArrayFlatten[Table[
    Which[i == j, sdiag[[i]], i == j + 1, slower[[j]],
      j == i + 1, ConjugateTranspose[slower[[i]]],
      True, zeroBlock[i, j]],
    {i, nBlocks}, {j, nBlocks}]];
  sigmaCSR = ArrayFlatten[Table[
    If[i == nBlocks && j == nBlocks, sigmaLast, zeroBlock[i, j]],
    {i, nBlocks}, {j, nBlocks}]];
  direct = Inverse[e scsr - hcsr - sigmaCSR];
  full = CentralGreenOverlap[e, {hcsr, scsr}, {sigmaCSR}];
  record["Complete CSR resolvent", Norm[full - direct, "Frobenius"]];

  ref = {{hdiag, hlower}, {sdiag, slower}};
  diagonal = CentralDiagonalBlockGreensOverlap[e, ref, {sigmaLast}];
  expectedDiagonal = Table[
    direct[[boundaries[[j]] + 1 ;; boundaries[[j + 1]],
      boundaries[[j]] + 1 ;; boundaries[[j + 1]]]],
    {j, nBlocks}];
  record["All CSR diagonal blocks", Max[MapThread[
    Norm[#1 - #2, "Frobenius"] &, {diagonal, expectedDiagonal}]]];

  lastBlock = CentralBlockGreensOverlap[e, ref, {sigmaLast}, "T"];
  record["Final CSR contact block", Norm[lastBlock - Last[expectedDiagonal],
    "Frobenius"]];
  lastColumn = Reverse[CentralBlockGreensOverlap[e, ref,
    {sigmaLast}, "LDOS"]];
  expectedColumn = Table[
    direct[[boundaries[[j]] + 1 ;; boundaries[[j + 1]],
      boundaries[[-2]] + 1 ;; boundaries[[-1]]]],
    {j, nBlocks}];
  record["Entire final CSR block column", Max[MapThread[
    Norm[#1 - #2, "Frobenius"] &, {lastColumn, expectedColumn}]]];

  (* 1b. Backward compatibility: orthogonal CSR input forms. *)
  orthogonalDirect = Inverse[
    e IdentityMatrix[nSites] - hcsr - sigmaCSR];
  orthogonalDiagonal = CentralDiagonalBlockGreens[
    e, {hdiag, hlower}, {sigmaLast}];
  record["Orthogonal CSR diagonal-block compatibility",
    Max[Table[Norm[orthogonalDiagonal[[j]] -
      orthogonalDirect[[boundaries[[j]] + 1 ;; boundaries[[j + 1]],
        boundaries[[j]] + 1 ;; boundaries[[j + 1]]]],
      "Frobenius"], {j, nBlocks}]]];
  orthogonalColumn = Reverse[CentralBlockGreens[
    e, {hdiag, hlower}, {sigmaLast}, "LCDV"]];
  record["Orthogonal CSR block-column compatibility",
    Max[Table[Norm[orthogonalColumn[[j]] -
      orthogonalDirect[[boundaries[[j]] + 1 ;; boundaries[[j + 1]],
        boundaries[[-2]] + 1 ;; boundaries[[-1]]]],
      "Frobenius"], {j, nBlocks}]]];
  record["Orthogonal CSR terminal compatibility",
    Norm[CentralBlockGreens[e, {hdiag, hlower},
      {sigmaLast}, "T"] - Last[orthogonalDiagonal], "Frobenius"]];

  (* 2. Lead modes 1..5 and interface self-energy at nonzero Im[e]. *)
  leadEnergy = 0.31 + 0.09 I;
  leadH0 = 0.17 herm[2];
  leadS0 = IdentityMatrix[2] + 0.01 herm[2];
  leadH1 = 0.12 rand[2, 2];
  leadS1 = 0.012 rand[2, 2];
  couplingH = 0.13 rand[2, 2];
  couplingS = 0.011 rand[2, 2];

  finiteLead = ArrayFlatten[Table[
    Which[
      i == j, leadEnergy leadS0 - leadH0,
      j == i + 1, leadEnergy leadS1 - leadH1,
      i == j + 1, leadEnergy ConjugateTranspose[leadS1] -
        ConjugateTranspose[leadH1],
      True, ConstantArray[0., {2, 2}]],
    {i, leadNumCells}, {j, leadNumCells}]];
  gdirect = Inverse[finiteLead][[1 ;; 2, 1 ;; 2]];

  Do[
    leadG = SurfaceGreenOverlap[leadEnergy,
      {{leadH0, leadH1}, {leadS0, leadS1}}, k];
    record["SurfaceGreenOverlap mode " <> ToString[k],
      Norm[leadG - gdirect, "Frobenius"], 1.*^-7],
    {k, 1, 5}];

  (* Orthogonal lead modes must retain their original calling convention. *)
  orthogonalLead = SurfaceGreen[leadEnergy, {leadH0, leadH1}, 1];
  Do[
    record["SurfaceGreen orthogonal mode " <> ToString[k],
      Norm[SurfaceGreen[leadEnergy, {leadH0, leadH1}, k] -
        orthogonalLead, "Frobenius"], 1.*^-7],
    {k, 2, 5}];
  record["Sigma orthogonal interface compatibility",
    Norm[Sigma[leadEnergy, {leadH0, leadH1, couplingH}, 1] -
      couplingH . orthogonalLead . ConjugateTranspose[couplingH],
      "Frobenius"], 1.*^-8];

  leadSigma = SigmaOverlap[leadEnergy,
    {{leadH0, leadH1, couplingH},
      {leadS0, leadS1, couplingS}}, 1];
  sigmaReference = (couplingH - leadEnergy couplingS) .
    gdirect . (ConjugateTranspose[couplingH] -
      leadEnergy ConjugateTranspose[couplingS]);
  record["SigmaOverlap mode 1", Norm[leadSigma - sigmaReference,
    "Frobenius"], 1.*^-7];

  Print["LGFF complex-energy regression results:"];
  Print[Dataset[report]];
  Print[If[AllTrue[Values[report], TrueQ[Lookup[#, "Pass"]] &],
    "ALL TESTS PASSED", "TEST FAILURE: inspect report"]];
]
