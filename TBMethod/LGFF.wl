(* ::Package:: *)

BeginPackage["TBMethod`LGFF`"]

(* Declare your package's public symbols here. *)

(* Exported symbols added here with SymbolName::usage *)
SurfaceGreen::usage = "surface Green function";
Sigma::usage = "Selfenergy from the one certain lead, source or drain.";
CentralGreen::usage = "Green's function for the central scattering area.";
Transmission::usage = "Transmission from Landauer-B\[UDoubleDot]ttiker formula.";
CentralBlockGreens::usage = "The blocks of the full Green's function in the last block column, for the (adaptively) partitioned central scattering region.";
LocalDOSRealSpace::usage = "Real-space local density of states from the Layered method.";
LocalDOSReciprocalSpace::usage = "Reciprocal-space local density of states.";
LocalCDV::usage = "Local current density vector field from the Layered method.";
TransportCoefficient::usage = "Calculates transmission or reflection coefficent by scattering matrix with Green's function method.";
CentralDiagonalBlockGreens::usage = "The diagonal blocks of the Green's function in the full form, corresponding to each layer from the (adaptively) partitioned central scattering region.";
HallAndLongitudinalResistances::usage = "Calculate the transverse Hall and longitudinal resistances simultaneously.";
HallAndLongitudinalConductances::usage = "Calculate the transverse Hall and longitudinal conductances simultaneously.";


Begin["`Private`"]
(* Implementation of the package *)
(*$DistributedContexts = {"Global`", "LGFF`"};*)

(* The transfer matrix *)
(*Tmat[matlist_] := Module[{positions(*, n = Length[matlist]*)},
	positions[m_] := Normal @ Partition[SparseArray[{{1} -> 1, {i_?OddQ} :> (i + 1)/2, {-1} -> 1}, 2 m, 2], 2];
	SparseArray[ Total[ Dot @@ Extract[matlist, positions[#] ] & /@ Range[Length[matlist]] ] ]
];*)

(* Another approach to find the inverse of a matrix instead of using `Inverse` directly *)
iden = IdentityMatrix[Length[#], SparseArray, WorkingPrecision -> MachinePrecision] &;
inverse[mat_, opts:OptionsPattern[LinearSolve]] := LinearSolve[mat, iden[mat], opts];

invbyQR[matrix_] :=
Module[{q, r},
	{q, r} = QRDecomposition[matrix];
	inverse[2][r] . q
];

(*transport_with_overlapping.nb*)
(*CSR as a whole*)
nonOrthoToOrthoTransf[e_][{HC_, SC_}] := HC + e(iden[HC] - SC);
(*surface Green's function of a lead*)
nonOrthoToOrthoTransf[e_][{{h0_, h1_}, {s0_, s1_}}] := {h0 + e(iden[h0] - s0), h1 - e s1};
(*self energy of a lead*)
nonOrthoToOrthoTransf[e_][{{h0_, h1_, V_}, {s0_, s1_, S_}}] := {h0 + e(iden[h0] - s0), h1 - e s1, V - e S};
(*CSR partitioned*)
nonOrthoToOrthoTransf[e_][{{ds_, os_}, {s0s_, s1s_}}] /; Length[ds] >= 2 := {ds + e(iden /@ ds - s0s), os - e s1s};

(* Surface Green function and selfenergy *)
(*SurfaceGreen[epsilon_, {h0_, h1_}] :=
Module[{id = iden[h0], inv = inverse[#, Method -> "Banded"] &, g0inverse, g0, t0, tttilde, tlist},
	g0inverse = epsilon id - h0;
	g0 = SparseArray[ inv[g0inverse] ];
	
	t0 = {g0.ConjugateTranspose[h1], g0.h1};
	tttilde[{a_, b_}] := Module[{(*tau,*) tauinversed},
		(*tau = id - (a.b + b.a);*)
		tauinversed = inv[ id - (a.b + b.a) ];
		{tauinversed.a.a, tauinversed.b.b}
		];
	tlist = FixedPointList[tttilde, t0, 2000];
	inv[ g0inverse - h1.Tmat[tlist] ]
];*)
(*SurfaceGreen[epsilon_, {h0_, h1_}, mode:(1|2):1] :=
Module[{id = iden[h0], inv = inverse[#, Method -> "Banded"] &, g0inverse, g0, t0, tttilde, M, S1, S2, n},
	g0inverse = epsilon id - h0;
	g0 = SparseArray[ inv[g0inverse] ];
	t0 = {g0 . ConjugateTranspose[h1], g0 . h1};
	
	T = Which[
		mode == 1,
		(tttilde[{a_, b_}] := Module[{tauinversed},
			tauinversed = inv[ id - (a . b + b . a) ];
			{tauinversed . a . a, tauinversed . b . b}
		];
		Tmat[ FixedPointList[tttilde, t0, 2000] ]),
		mode == 2,
		(M = ArrayFlatten[({{#, -# . First[t0]}, {id, 0.}})& [inverse[Last @ t0]]] // SparseArray; n = Length[id];
		(*{S1, S2} = Partition[SortBy[Eigensystem[M, Method \[Rule] "Direct"]\[Transpose], Abs@*First]\[LeftDoubleBracket];;n, 2\[RightDoubleBracket]\[Transpose], n];*)
		{S1, S2} = Partition[Eigenvectors[M, -n, Method -> "Direct"]\[Transpose], n];
		S1 . inverse[S2])
	];
	
	inv[ g0inverse - h1 . T ]
];*)
(*SurfaceGreen[epsilon_, {h0_, h1_}, mode:(1|2):1] :=
Module[{id = iden[h0], inv = inverse[#, Method -> "Banded"] &, g0inverse, g0, t0, tttilde, TTtilde, T, M, S1, S2, n},
	g0inverse = epsilon id - h0;
	
	T = Which[
		(*iterative 2^n method*)
		mode == 1,
		g0 = SparseArray[ inv[g0inverse] ];
		t0 = {g0 . h1\[ConjugateTranspose], g0 . h1};
		(tttilde[{a_, b_}] := Module[{tauinversed},
			tauinversed = inv[ id - (a . b + b . a) ];
			{tauinversed . a . a, tauinversed . b . b}
		];
		TTtilde[{{t_, tt_}, {T_, Tt_}}] := Module[{newt, newtt},
			{newt, newtt} = tttilde[{t, tt}];
			{{newt, newtt}, {T + Tt . newt, Tt . newtt}}
		];
		FixedPoint[TTtilde, {t0, t0}, 2000][[2, 1]]
		),
		(*transfer matrix method*)
		mode == 2,
		((*M = ArrayFlatten[({{#, -# . First[t0]}, {id, 0.}})& [inverse[Last @ t0]]] // SparseArray;*)
		M = ArrayFlatten[({{# . g0inverse, -# . h1\[ConjugateTranspose]}, {id, 0.}})& [inverse[h1]]] // SparseArray;
		n = Length[id];
		(*{S1, S2} = Partition[SortBy[Eigensystem[M, Method \[Rule] "Direct"]\[Transpose], Abs@*First]\[LeftDoubleBracket];;n, 2\[RightDoubleBracket]\[Transpose], n];*)
		{S1, S2} = Partition[Eigenvectors[M, -n, Method -> "Direct"]\[Transpose], n];
		S1 . inverse[S2])
	];
	
	inv[ g0inverse - h1 . T ]
];*)


(*SurfaceGreenOverlap[e_, {{h0_, h1_}, {s0_, s1_}}, mode:(1|2|3):1] := SurfaceGreen[e, nonOrthoToOrthoTransf[e][{{h0, h1}, {s0, s1}}], mode];*)
SurfaceGreen[e_, {h0_, h1_}, {s0_, s1_}, mode:(1|2|3):3] := SurfaceGreen[e, nonOrthoToOrthoTransf[e][{{h0, h1}, {s0, s1}}], mode];
SurfaceGreen[e_, {h0_, h1_}, mode:(1|2|3):3] :=
Module[{id = iden[h0], inv = inverse[#, Method -> "Banded"] &, g0inverse, g0, t0, tttilde, TTtilde, T, MH, S1, S2, n},
	g0inverse = e id - h0;
	
	T = Which[
		(*iterative 2^n method*)
		mode == 1,
		g0 = SparseArray[ inv[g0inverse] ];
		t0 = {g0 . h1\[ConjugateTranspose], g0 . h1};
		(tttilde[{a_, b_}] := Module[{tauinversed},
			tauinversed = inv[ id - (a . b + b . a) ];
			{tauinversed . a . a, tauinversed . b . b}
		];
		TTtilde[{{t_, tt_}, {T_, Tt_}}] := Module[{newt, newtt},
			{newt, newtt} = tttilde[{t, tt}];
			{{newt, newtt}, {T + Tt . newt, Tt . newtt}}
		];
		FixedPoint[TTtilde, {t0, t0}, 2000][[2, 1]]
		),
		(*transfer matrix method*)
		mode == 2 || mode == 3,
		n = Length[id];
		(MH = If[mode == 2,
			SparseArray @ ArrayFlatten[({{# . g0inverse, -# . h1\[ConjugateTranspose]}, {id, 0.}}) & [inv[h1]]],
			SparseArray @* ArrayFlatten /@ {{{g0inverse, -h1\[ConjugateTranspose]}, {id, 0.}}, {{h1, 0.}, {0., id}}}
		];
		{S1, S2} = Partition[SortBy[Eigensystem[MH, Method -> "Direct"]\[Transpose], Abs @* First][[;;n, 2]]\[Transpose], n];
		(*{S1, S2} = Partition[Eigenvectors[MH, -n, Method -> "Direct"]\[Transpose], n];*)
		S1 . inv[S2])
	];
	
	inv[ g0inverse - h1 . T ]
];


(*Sigma[epsilon_, {h0_, h1_, H01_}, mode:(1|2|3):1] :=
Module[{gsurface},
	gsurface = SurfaceGreen[epsilon, {h0, h1}, mode];
	SparseArray[H01 . gsurface . H01\[ConjugateTranspose]]
];*)

(*SigmaOverlap[e_, {{h0_, h1_, H01_}, {s0_, s1_, S01_}}, mode:(1|2|3):1] := Sigma[e, nonOrthoToOrthoTransf[e][{{h0, h1, H01}, {s0, s1, S01}}], mode];*)
Sigma[e_, {h0_, h1_, H01_}, {s0_, s1_, S01_}, mode:(1|2|3):3] := Sigma[e, nonOrthoToOrthoTransf[e][{{h0, h1, H01}, {s0, s1, S01}}], mode];
Sigma[e_, {h0_, h1_, H01_}, mode:(1|2|3):3] :=
Module[{gsurface, len = Length[H01], \[Eta] = Im[e]},
	If[Norm[H01, "Frobenius"] >= \[Eta], 
		gsurface = SurfaceGreen[e, {h0, h1}, mode]; SparseArray[H01 . gsurface . H01\[ConjugateTranspose]],
		SparseArray[{}, {1, 1} len]
	]
];

(*Green's function of the central region*)
(*CentralGreenOverlap[e_,  {HC_, SC_}, Sigmas_List, method_:1] := CentralGreen[e, nonOrthoToOrthoTransf[e][{HC, SC}], Sigmas, method];*)
CentralGreen[e_,  HC_, SC_, Sigmas_List, method_:1] := CentralGreen[e, nonOrthoToOrthoTransf[e][{HC, SC}], Sigmas, method];
CentralGreen[e_, HC_, Sigmas_List, method_:1] :=
Module[{id = iden[HC], GCinverse, inv = inverse[#, Method -> "Multifrontal"] &},
	GCinverse = e id - HC;
	(*SparseArray[ #[GCinverse - Total @ Sigmas] ] & @ Which[method == 1, inverse[#, Method -> "Multifrontal"]&, method == 2, invbyQR]*)
	Which[method == 1, inv, method == 2, invbyQR][GCinverse - Total @ Sigmas]
];

(*CentralBlockGreensOverlap[e_, {blockHs:{ds_, os_}, blockSs:{s0s_, s1s_}}, sigmas_, mode:("T"|"LDOS"|"LCDV"):"T"] := CentralBlockGreens[e, nonOrthoToOrthoTransf[e][{blockHs, blockSs}], sigmas, mode];*)
CentralBlockGreens[e_, blockHs:{ds_, os_}, blockSs:{s0s_, s1s_}, sigmas_, mode:("T"|"LDOS"|"LCDV"):"T"] := CentralBlockGreens[e, nonOrthoToOrthoTransf[e][{blockHs, blockSs}], sigmas, mode];	
CentralBlockGreens[e_, blockHs:{ds_, os_}, sigmas_, mode:("T"|"LDOS"|"LCDV"):"T"] /; (Subtract @@ (Length /@ blockHs) == 1) := 
Module[{inv, diagblocks, iteratefunc, arguments},
	inv = inverse[#, Method -> "Multifrontal"] &;
	diagblocks = e iden[#] - # & /@ MapAt[Total[sigmas] + # &, ds, -1];
	arguments = Sequence[iteratefunc, inv[First[diagblocks]], Transpose[{-os, Rest @ diagblocks}] ];
	iteratefunc = inv[ Last[#2] - First[#2] . # . First[#2]\[ConjugateTranspose] ] &;
	Which[
		mode == "T", Fold[arguments],
		mode == "LDOS" || mode == "LCDV",
		Module[{Fi, Fioidag, iteratefunc2},
			Fi = FoldList[arguments];
			Fioidag = MapThread[Dot, {Most[Fi], -ConjugateTranspose /@ os}];
			iteratefunc2 = -#2 . # &;
			FoldList[iteratefunc2, Last[Fi], Reverse[Fioidag]]
		],
		True, 0
	]
];

(*The diagonal blocks of the Green's functions*)
(*CentralDiagonalBlockGreensOverlap[e_, {blockHs:{ds_, os_}, blockSs:{s0s_, s1s_}}, sigmas_] := CentralDiagonalBlockGreens[e, nonOrthoToOrthoTransf[e][{blockHs, blockSs}], sigmas];*)
CentralDiagonalBlockGreens[e_, blockHs:{ds_, os_}, blockSs:{s0s_, s1s_}, sigmas_] := CentralDiagonalBlockGreens[e, nonOrthoToOrthoTransf[e][{blockHs, blockSs}], sigmas];
CentralDiagonalBlockGreens[e_, blockHs:{ds_, os_}, sigmas_] /; And[Length[ds] >= 2, Subtract @@ (Length /@ blockHs) == 1] :=
Module[{inv, diagblocks, blockGinvs, blockGinvsrev, iteratefunc, Fisinner, Fisouter, foldlist},
	inv = inverse[#1, Method -> "Multifrontal"] &;
	foldlist = FoldList[Last[#2] - First[#2] . inv[#] . First[#2]\[ConjugateTranspose] &];
	diagblocks = e iden[#] - # & /@ MapAt[Total[sigmas] + # &, ds, -1];
	
	blockGinvs = {diagblocks, -os};
	Fisinner = foldlist[blockGinvs[[1, 1]], {blockGinvs[[2]], Rest @ blockGinvs[[1]]}\[Transpose]];
	
	blockGinvsrev = Reverse[{diagblocks, -ConjugateTranspose /@ os(*ConjugateTranspose[os, 2 <-> 3]*)}, 2];
	Fisouter = foldlist[blockGinvsrev[[1, 1]], {blockGinvsrev[[2]], Rest @ blockGinvsrev[[1]]}\[Transpose]] // Reverse;
	
	MapThread[inv[# - #2] &][{Fisinner + Fisouter, diagblocks}]
];



Transmission[GCSR_, \[CapitalSigma]s_] :=
Module[{\[CapitalGamma]p, \[CapitalGamma]q},
	{\[CapitalGamma]p, \[CapitalGamma]q} = I (# - #\[ConjugateTranspose]) & /@ \[CapitalSigma]s;
	(*{\[CapitalGamma]p, \[CapitalGamma]q} = -2 Im /@ \[CapitalSigma]s;*)(*This is very wrong!*)
	Tr[ \[CapitalGamma]p . GCSR . \[CapitalGamma]q . GCSR\[ConjugateTranspose] ]
];

scatteringMatrix[GCSR_, \[CapitalSigma]s_, zero_ : 1.*^-4] :=
Module[{\[CapitalGamma], \[CapitalGamma]sqrt, \[CapitalGamma]s, \[CapitalGamma]ssqrt, len = Length[GCSR], linewidth = I (# - #\[ConjugateTranspose]) &, sqrtmat = MatrixPower[#, 1/2] &, Mtot},
	If[Equal @@ \[CapitalSigma]s,
		(\[CapitalGamma] = linewidth[\[CapitalSigma]s[[2]]]; \[CapitalGamma]sqrt = sqrtmat[\[CapitalGamma]];
		Mtot = MatrixRank[\[CapitalGamma]\[ConjugateTranspose] . \[CapitalGamma], Tolerance -> zero, ZeroTest -> (Abs[#]^2 < zero &)];
		{- IdentityMatrix[len] + I \[CapitalGamma]sqrt . GCSR . \[CapitalGamma]sqrt, - len + Mtot}),
		(\[CapitalGamma]s =  linewidth /@ \[CapitalSigma]s; \[CapitalGamma]ssqrt = sqrtmat /@ \[CapitalGamma]s;
		{I # . GCSR . #2 & @@ \[CapitalGamma]ssqrt, 0})
	]
];

TransportCoefficient[GCSR_, \[CapitalSigma]s_ , zero_ : 1.*^-4] :=
Module[{matSq = #\[ConjugateTranspose] . # &, smat, Mtot},
	{smat, Mtot} = scatteringMatrix[GCSR, \[CapitalSigma]s, zero];
	Tr[matSq[smat]] + Mtot
];

(*LocalDOSRealSpace[Gs_, sigmas_, layeredpts_Association, innerdof_:1] :=
Module[{gamma, innerdofldos, ldos},
	gamma = I (# - ConjugateTranspose[#]) & @ Total[sigmas];
	(*gamma = Im @ Total[sigmas];*)
	innerdofldos = (*-*)1/\[Pi] Diagonal[ # . gamma . ConjugateTranspose[#] ] & /@ Gs;
	(*innerdofldos = -1/\[Pi] Diagonal[ Im[#] ] & /@ Gs;*) (*WRONG!*)
	ldos = BlockMap[Total, #, innerdof] & /@ innerdofldos;
	MapThread[Append, Join @@@ {Values[layeredpts], Reverse[ldos]}]
];*)
LocalDOSRealSpace[innerproj_, layeredpts_Association][Gs_, sigmas_] :=
Module[{gamma, innerdofldos, ldos, \[Rho]s, innerdof = Dimensions[innerproj]},
	gamma = I (# - #\[ConjugateTranspose]) & @ Total[sigmas];
	(*gamma = Im @ Total[sigmas];*) (* WRONG!*)
	innerdofldos = 1/(2\[Pi]) Diagonal[Partition[# . gamma . #\[ConjugateTranspose], innerdof]] & /@ Gs;
	ldos = Map[Tr[innerproj . #] &, innerdofldos, {2}];
	MapThread[Append, Join @@@ {Values[layeredpts], Reverse[ldos]}]
];

(*LocalDOSReciprocalSpace[{k_, \[Epsilon]_}, {h00_, h01_}, mode:(1|2):1] := LocalDOSReciprocalSpace[{k, \[Epsilon]}, {h00, h01, h00}, mode];
LocalDOSReciprocalSpace[{k_, \[Epsilon]_}, {h00_, h01_, H00_}, mode:(1|2):1] :=
Module[{zero = 1.*^-4, \[CapitalSigma]},
	\[CapitalSigma] = Sigma[Complex[\[Epsilon], zero], {h00, h01, h01}];
	-Im @ Tr @ CentralGreen[Complex[\[Epsilon], zero], H00, {\[CapitalSigma]}]
];*)
(*LocalDOSReciprocalSpace[\[Epsilon]_, {HLeadBloch_, HLead12_}, mode:(1|2|3):1] := LocalDOSReciprocalSpace[\[Epsilon], {HLeadBloch, HLead12}, {HLeadBloch, HLead12}, mode];
LocalDOSReciprocalSpace[\[Epsilon]_, {HLeadBloch_, HLead12_}, HCSRBloch_, mode:(1|2|3):1] := LocalDOSReciprocalSpace[\[Epsilon], {HLeadBloch, HLead12}, {HCSRBloch, HLead12}, mode];
LocalDOSReciprocalSpace[\[Epsilon]_, {HLeadBloch_, HLead12_}, {HCSRBloch_, HCSRLead1_}, mode:(1|2|3):1] :=
Module[{zero = 1.*^-4, \[CapitalSigma]},
	\[CapitalSigma] = Sigma[Complex[\[Epsilon], zero], {HLeadBloch, HLead12, HCSRLead1}, mode];
	-Im @ Tr @ CentralGreen[Complex[\[Epsilon], zero], HCSRBloch, {\[CapitalSigma]}]
];*)
LocalDOSReciprocalSpace[\[Epsilon]_, {HLeadBloch_?MatrixQ, HLead12_?MatrixQ}, mode:(1|2|3):3] := LocalDOSReciprocalSpace[\[Epsilon], {{HLeadBloch, HLead12}, {HLeadBloch, HLead12}}, mode];
LocalDOSReciprocalSpace[\[Epsilon]_, {HLeadBloch_?MatrixQ, HLead12_?MatrixQ}, HCSRBloch_?MatrixQ, mode:(1|2|3):3] := LocalDOSReciprocalSpace[\[Epsilon], {{HLeadBloch, HLead12}, {HCSRBloch, HLead12}}, mode];
LocalDOSReciprocalSpace[\[Epsilon]_, {{HLeadBloch_?MatrixQ, HLead12_?MatrixQ}, {HCSRBloch_?MatrixQ, HCSRLead1_?MatrixQ}}, mode:(1|2|3):3] :=
Module[{zero = 1.*^-4, \[CapitalSigma]},
	\[CapitalSigma] = Sigma[Complex[\[Epsilon], zero], {HLeadBloch, HLead12, HCSRLead1}, mode];
	-Im @ Tr @ CentralGreen[Complex[\[Epsilon], zero], HCSRBloch, {\[CapitalSigma]}]
];


(*Local current density vector*)
(*two auxilliary functions*)
currentTensorBlocks[Gs_, sigmas_, blockHs: {ds_, os_}, innerdof_:1] :=
Module[{gamma, blockGns, Gsre = Reverse[Gs], jblock0, jblock0innersummed},
	gamma = I (# - #\[ConjugateTranspose]) & [Total[sigmas]];
	blockGns = {# . gamma . #\[ConjugateTranspose] & /@ Gsre, # . gamma . #2\[ConjugateTranspose] & @@@ Partition[Gsre, 2, 1]};
	(*jblock0 = -Im[MapAt[ConjugateTranspose, {2, All}][blockHs] blockGns];*)
	(*jblock0 = Im[MapAt[Transpose, {2, All}][blockHs] blockGns];*)
	jblock0 = Im[Map[Transpose, blockHs, {2}] blockGns];(*!!!*)
	jblock0innersummed = Table[BlockMap[Total[#, 2] &, #, {1, 1}innerdof] & /@ x, {x, jblock0}];
	(*How to sum the internal degree of freedom?*)
	Append[jblock0innersummed, -Transpose /@ jblock0innersummed[[2]]]
];(*bond current in layered block form*)

jvecfield[ptsCSR_, jmat_] := jvecfield[{ptsCSR, ptsCSR}, jmat];
jvecfield[{ptsf_, ptsi_}, jmat_] := Module[{jvec, groupfunc, zero = 1.*^-3},
	jvec[find_, iindslengths_] := Module[{ptf = ptsf[[find]], vec},
		(*vec = Sum[Normalize[ptf - ptsi[[Keys[il]]]]Values[il], {il, iindslengths}]; (*Normalize[] might be unnecessary.*)*)
		(*vec = Sum[(*Normalize@*)(ptf - ptsi[[Keys[il]]]) Values[il], {il, iindslengths}];*)
		vec = Sum[(*Normalize@*)(ptsi[[Keys[il]]] - ptf) Values[il], {il, iindslengths}];
		(*If[Norm[vec] < zero, Nothing, {ptf, vec}]*)
		(*If[Norm[vec] < zero, Nothing, ptf -> vec]*)
		ptf -> vec
	];
	groupfunc = First @* Keys -> (Last @ Keys[#] -> Values[#] &);
	KeyValueMap[jvec] @ GroupBy[Most @ ArrayRules[jmat], groupfunc]
];(*from bond current to current field*)

jtensorFromBlocks[jtensorblocks:{a_, _, _}] /; (Equal @@ (Length /@ jtensorblocks - {1, 0, 0})) :=
Module[{inds, inddiag, indupoffd, n = Length[a], blocks = Join @@ jtensorblocks},
	inddiag = Table[{i, i}, {i, n}]; indupoffd = Table[{i, i+1}, {i, n-1}];
	inds = Join[inddiag, indupoffd, Reverse[indupoffd, 2]];
	SparseArray`SparseBlockMatrix[Thread[inds -> blocks]]
];(*from layered blocks to a whole block*)

(*LocalCDV[ptslayered_Association, currenttensorblocks_, innerdof_:1] :=
Module[{ptscsr = Values[ptslayered], innersummed, ptspairs, ptspairsfinal, js},
	ptspairs = Partition[ptscsr, 2, 1]; ptspairsfinal = {ptscsr, ptspairs, Reverse[ptspairs, 2]};
	innersummed = Table[BlockMap[Total[#, 2] &, #, {1, 1}innerdof] & /@ x, {x, currenttensorblocks}];
	js = Join @@ Table[Join @@ MapThread[jvecfield, {ptspairsfinal[[i]], innersummed[[i]]}], {i, 3}];
	KeyValueMap[List] @ (Total /@ GroupBy[js, First -> Last])
];*)

(*LocalCDV[Gs_, sigmas_, blockHs: {ds_, os_}, layeredpts_Association, innerdof_:1] :=
Module[{ptscsr = Values[layeredpts], innersummed, ptspairs, ptspairsfinal, js, currenttensorblocks},
	currenttensorblocks = currentTensorBlocks[Gs, sigmas, layeredpts, blockHs];
	ptspairs = Partition[ptscsr, 2, 1]; ptspairsfinal = {ptscsr, ptspairs, Reverse[ptspairs, 2]};
	innersummed = Table[BlockMap[Total[#, 2] &, #, {1, 1}innerdof] & /@ x, {x, currenttensorblocks}];
	js = Join @@ Table[Join @@ MapThread[jvecfield, {ptspairsfinal[[i]], innersummed[[i]]}], {i, 3}]; (*this is wrong*)
	KeyValueMap[List] @ (Total /@ GroupBy[js, First -> Last])
];*)

LocalCDV[Gs_, sigmas_, blockHs: {ds_, os_}, layeredpts_Association, innerdof_:1] :=
Module[{ptscsr = Join @@ layeredpts, innersummed, js, currenttensorblocks, jtensorfull},
	currenttensorblocks = currentTensorBlocks[Gs, sigmas, blockHs, innerdof];
	jtensorfull = jtensorFromBlocks[currenttensorblocks];
	js = jvecfield[ptscsr, jtensorfull];
	(*KeyValueMap[List] @ (Total /@ GroupBy[js, First -> Last])*)
	KeyValueMap[List] @* Merge[Total] @ js
];

(*For five and six T cases; The numbering of leads follows the order of from left to right and from up to down.
And the source and drain should be numbered at the beginning and the ending.*)
(*
     2   
1 ------- 4
     3

   2   4
1 ------- 5
   3

   2   4
1 ------- 6
   3   5
*)
HallAndLongitudinalResistances[\[Epsilon]_, hcsrdod_, leadshs_] :=
Module[{\[ScriptCapitalT], \[ScriptCapitalT]func, Rfunc, cnup = 1.*^7, ter = Length[leadshs], transmissions},
	\[ScriptCapitalT]func = # - DiagonalMatrix[Total[#, {2}]] &;
	Rfunc = # - {#2, #3} & @@ Rest[LinearSolve[#][UnitVector[ter - 1, 1]]] &;
	transmissions = Module[{\[CapitalSigma]s, blockG},
		\[CapitalSigma]s = Sigma[\[Epsilon], #, 3] & /@ leadshs;
		blockG = CentralBlockGreens[\[Epsilon], hcsrdod, \[CapitalSigma]s, "T"];
		Table[If[p == q || p == ter, 0., Transmission[blockG, \[CapitalSigma]s[[{p, q}]]]], {p, ter}, {q, ter}]
	];
	\[ScriptCapitalT] = Drop[\[ScriptCapitalT]func[-transmissions], -1, -1];
	If[LinearAlgebra`Private`MatrixConditionNumber[\[ScriptCapitalT]] > cnup, {"NaN", "NaN"},
		Rfunc[\[ScriptCapitalT]]
	]
];

HallAndLongitudinalConductances[\[Epsilon]_, hcsrdod_, leadshs_] :=
Module[{RH, RL},
	{RH, RL} = HallAndLongitudinalResistances[\[Epsilon], hcsrdod, leadshs];
	{RH, RL}/(RH^2 + RL^2)
];

(*HallAndLongitudinalConductances[\[Epsilon]_, hcsrdod_, leadshs_] :=
Module[{\[ScriptCapitalT], comat, \[ScriptCapitalT]func, Rfunc, RH, RL, inverse, cnup = 1.*^7, ter = Length[leadshs], transmissions},
	\[ScriptCapitalT]func = # - DiagonalMatrix[Total[#, {2}]] &;
	Rfunc = # - {#2, #3} & @@ Rest[LinearSolve[#][UnitVector[ter - 1, 1]]] &;
	transmissions = Module[{\[CapitalSigma]s, blockG},
		\[CapitalSigma]s = Sigma[\[Epsilon], #, 3] & /@ leadshs;
		blockG = CentralBlockGreens[\[Epsilon], hcsrdod, \[CapitalSigma]s, "T"];
		Table[If[p == q || p == ter, 0., Transmission[blockG, \[CapitalSigma]s[[{p, q}]]]], {p, ter}, {q, ter}]
	];
	\[ScriptCapitalT] = Drop[\[ScriptCapitalT]func[-transmissions], -1, -1];
	If[LinearAlgebra`Private`MatrixConditionNumber[\[ScriptCapitalT]] > cnup, {"NaN", "NaN"},
		{RH, RL} = Rfunc[\[ScriptCapitalT]]; {RH, RL}/(RH^2 + RL^2)
	]
];*)


End[] (* End `Private` *)

EndPackage[]
