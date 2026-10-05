(* ::Package:: *)

(* ::Text:: *)
(*This file collects old functions as an archive, NOT for the purpose to run.*)
(**)
(*A function being old here almost means it is created earlier in the D&R phase but works well, however not as efficient as the latest version in the runtime files.*)
(**)
(*But hidden bugs and/or errors cannot be 100% ruled out. And because they are old, I will not devote more effort to repair them.*)
(**)
(*Or, at least the newer versions are their efficiency reparations.*)


(* ::Title::Closed:: *)
(*MDConstruct`**)


(*PhaseFactor2DAB[B_, \[Phi]A_][ptf:{_, _, _}, pti:{_, _, _}] := PhaseFactor2DAB[B, \[Phi]A][Most @ ptf, Most @ pti]
PhaseFactor2DAB[B_, \[Phi]A_][ptf:{_, _}, pti:{_, _}] :=
Module[{xi, yi, xj, yj, \[CurlyPhi]},
	{{xi, yi}, {xj, yj}} = {ptf, pti};
	(*B \[Pi] ( xj yi - xi yj - (xi yi - xj yj) Cos[2\[Phi]A] + (xi^2 - xj^2 - yi^2 + yj^2) Sin[2\[Phi]A]/2)*)
	(*electron has a negative charge -e*)
	\[CurlyPhi] = B \[Pi] (- xj yi + xi yj + (xi yi - xj yj) Cos[2\[Phi]A] + (-xi^2 + xj^2 + yi^2 - yj^2) Sin[2\[Phi]A]/2);
	Exp[I \[CurlyPhi]]
];*)

(*coordspattern = {{__?NumericQ}..|{{__?NumericQ}, True|False}..|(_ -> {__?NumericQ})..};*)
(*coordspattern = {{__?NumericQ}..|(_ -> {__?NumericQ})..};*)
(*coordspattern={{__?NumericQ}..}|{Rule[_, {__?NumericQ}]..};*)

(*neighbourInfos[pts:coordspattern, dist_Real] := neighbourInfos[{pts, pts}, dist];*)(*This kind of pattern splitting leads to ambiguity.*)

(*HMatrixFromHoppings[pts:coordspattern, tFunc_, dist_Real, sumfunc_:Sum] := HMatrixFromHoppings[{pts, pts}, tFunc, dist, sumfunc];*)
(*General::maxrec: "Recursion limit exceeded; positive match might be missed."*)

(*HMatrixFromHoppings[fipts:{fpts:coordspattern, ipts:coordspattern}, tFunc_, dist_Real, sumfunc_:Sum] :=
Block[{dim = Length /@ fipts, neighbourinfos = neighbourInfos[fipts, dist], len, (*lenup = 1000,*) summand},
	len = Length[neighbourinfos];
	summand[ij_] := KroneckerProduct[SparseArray[ij -> 1, dim], tFunc[fpts[[#]], ipts[[#2]]] & @@ ij];
	(*Which[
		len == 0, KroneckerProduct[SparseArray[{}, dim], tFunc[fpts[[1]], ipts[[1]]]],
		lenup > len > 0, Sum[summand[ij], {ij, neighbourinfos}],
		lenup <= len, ParallelSum[summand[ij], {ij, neighbourinfos}]
	]*)
	Which[
		len == 0, KroneckerProduct[SparseArray[{}, dim], tFunc[fpts[[1]], ipts[[1]]]],
		len > 0 , sumfunc[summand[ij], {ij, neighbourinfos}]
	]
];*)
(*HMatrixFromHoppings[pts:coordspattern, tFunc_, dist_Real] := HMatrixFromHoppings[{pts, pts}, tFunc, dist];
HMatrixFromHoppings[fipts:{fpts:coordspattern, ipts:coordspattern}, tFunc_, dist_Real] :=
Module[{dim = Length /@ fipts, fptsnobool, iptsnobool, nfunc, neighbourinfos, fillingrules, a0 = 1},
	{fptsnobool, iptsnobool} = If[FreeQ[#, True|False], #, #[[;;, 1]]]& /@ fipts;
	(*nfunc is a NearestFunction with respect to the group of initial points, which should constitute the column index*)
	nfunc = Nearest[iptsnobool -> Automatic, WorkingPrecision -> MachinePrecision(*Method -> Automatic(*{"KDtree","LeafSize"->50}*)*)];
	(*Within the disk/ball with a certain radius (= dist a0) centered at each final point, find the encompassed initial points; each point is given as its index in its corresponding ordered set (final point in fpts & initial point in ipts).*)
	neighbourinfos = Flatten[MapIndexed[Thread[{#2[[1]], #}]&, nfunc[fptsnobool, {All, dist a0}]], 1];(*{{j..}..} -> {{i, {j..}}..} -> {{i, j}..}*)
	
	(*Summation saves RAM considerably. Fine tunability lies in the methods of Parallelize[]*)
	If[Length[neighbourinfos] < 1000, Sum, ParallelSum][
		KroneckerProduct[SparseArray[ij -> 1, dim], tFunc[fpts[[#]], ipts[[#2]]]& @@ ij],
		{ij, neighbourinfos}
	]
];*)

(*HBloch[vk_, h0010s:<|({__?NumericQ} -> _SparseArray)..|>] :=
Module[{hermitize = # + #\[HermitianConjugate] &},
	First[h0010s] + hermitize[KeyValueMap[Exp[I # . vk] #2 &, Rest[h0010s]] // Total]
];

HBlochFull[vk_, vecaHa_Association] := Total[KeyValueMap[Exp[I # . vk] #2 &, vecaHa]];*)

(*HPDTensorsBloch[n_Integer?NonNegative][vk_, h0010s : <|({__?NumericQ} -> _SparseArray) ..|>] :=
Module[{fullassoc},
	fullassoc = Join[h0010s, Association @ KeyValueMap[(-#1 -> #2\[HermitianConjugate]) &, Rest[h0010s]]];
	HPDTensorsBlochFull[n][vk, fullassoc]
];*)

(*ParallelHCSRDiagOffDiagBlocks[CSRptsgrouped_Association, tFunc_, dup_] :=
Module[{hdfill, hofill, hdblocks, hoblocks, csrpts = Values[CSRptsgrouped], len = Length[CSRptsgrouped]},
	If[len == 1, HMatrixFromHoppings[Join[csrpts, csrpts], tFunc, dup],
		(hdfill = p |-> HMatrixFromHoppings[{p, p}, tFunc, dup];
		hofill = p |-> HMatrixFromHoppings[p, tFunc, dup];
		hdblocks = hdfill ~ParallelMap~ csrpts;
		hoblocks = hofill ~ParallelMap~ Transpose[{Rest[csrpts], Most[csrpts]}];
		(*hoblocks = Parallelize[MapThread[hofill @* List, {Rest[csrpts], Most[csrpts]}]];*)
		{hdblocks, hoblocks})
	]
];*)

(*HLeadBlocks[CSRptsgrouped_Association, tFunc_, dup_, leadpts_] := Table[HMatrixFromHoppings[fipts, tFunc, dup], {fipts, {{#[[1]], #[[1]]}, #, {Values[CSRptsgrouped][[-1]], #[[1]]}} & [leadpts]}];*)

(*Points in CSR without labels*)
(*AdaptivePartition[{ptslead1stcell_, ptscsr:coordspattern[[{1}, 1]]}, dup_] :=
Module[{groupfunc, iterate, ptscsrgrouped},
	groupfunc[{x_, y_}] := Module[{ptsneighbors, layerneighbors},
		ptsneighbors = Join @@ NearestTo[x, {All, dup}, WorkingPrecision -> MachinePrecision, Method -> "KDTree"][y];
		layerneighbors = DeleteDuplicates[ptsneighbors];
		{x, {#, Complement[y, #]}} & [layerneighbors]
	];
	iterate = FlattenAt[-1] @* SubsetMap[groupfunc, -2;;];
	(*ptscsrgrouped = NestWhile[iterate, {ptslead1stcell, ptscsr}, Last[#] != {} &][[2;;-2]];*)
	ptscsrgrouped = Rest @ NestWhile[iterate, {ptslead1stcell, ptscsr}, Last[#] != {} &, 1, \[Infinity], -1];
	<|MapIndexed[#2[[1]] -> # &, Reverse[ptscsrgrouped]]|>
];
(*Points in CSR with labels: this second version actually works for both patterns, but slower than the first version for the first pattern.*)
AdaptivePartition[ptsleadcsr:{ptslead1stcell_, ptscsr:coordspattern[[{1}, 2]]}, dup_] :=
Module[{neighborindex, indexneighbortolead, indexgrouped, iterate, ptsleadnobool, ptscsrnobool, indexlayered, len = Length[ptscsr], nf},
	{ptsleadnobool, ptscsrnobool} = If[FreeQ[#, Rule[_, _]], #, Values[#]] & /@ ptsleadcsr;
	nf = Nearest[ptscsrnobool -> Automatic, WorkingPrecision -> MachinePrecision, Method -> "KDTree"];
	neighborindex[x_] := DeleteDuplicates @* Join @@ nf[x, {All, dup}];
	indexneighbortolead = neighborindex[ptsleadnobool];
	iterate = neighborindex[ptscsrnobool[[#]]] &;
	indexgrouped = Reverse @ NestWhileList[iterate, indexneighbortolead, Length[#] < len &];
	indexlayered = BlockMap[Complement @@ # &, Append[indexgrouped, {}], 2, 1];
	<|MapIndexed[#2[[1]] -> ptscsr[[#]] &, indexlayered]|>
]*)

(*AdaptivePartition[ptsleadcsr:{ptslead1stcell:coordspattern, ptscsr:coordspattern}, dup_, opts: OptionsPattern[Nearest]] :=
Module[{neighborindex, nf, indexgrouped, iterate, ptsleadnobool, ptscsrnobool, indexlayered, initial},
	{ptsleadnobool, ptscsrnobool} = If[FreeQ[#, Rule[_, _]], #, Values[#]] & /@ ptsleadcsr;
	nf = Nearest[ptscsrnobool -> Automatic, opts, WorkingPrecision -> MachinePrecision, Method -> "KDTree"];
	neighborindex[x_] := DeleteDuplicates @* Join @@ nf[x, {All, dup}];
	iterate = Append[#, Complement[neighborindex[ptscsrnobool[[Last[#]]]], #[[-1]], #[[-2]]]] &;
	initial = {{}, neighborindex[ptsleadnobool]};
	indexlayered = Reverse @ Rest @ NestWhile[iterate, initial, Last[#] != {} &, 1, \[Infinity], -1];
	<|MapIndexed[#2[[1]] -> ptscsr[[#]] &, indexlayered]|>
];*)

(*HBlochsForSpecFunc[vk_, ptscellsvas_, tfunc_, dup_] :=
Module[{HLeadIntraInterReal, HCSRIntraLeadInterReal, fillfunc, keys = Keys[ptscellsvas], len = Length[ptscellsvas]},
	fillfunc[indfs_, indi_] := HMatrixFromHoppings[{#, ptscellsvas[[indi]]}, tfunc, dup] & /@ KeyMap[# - keys[[indi]] &][ptscellsvas[[indfs]]];
	(*fillfunc[indfs_, indi_] := HMatrixFromHoppings[{#, ptscellsvas[[indi]]}, tfunc, dup] & /@ ptscellsvas[[indfs]];*)
	HLeadIntraInterReal = fillfunc @@@ {{2;;4, 3}, {2;;4, 1}};
	Which[
		len == 4, HBlochFull[vk, #]& /@ HLeadIntraInterReal,
		len == 7,
		(HCSRIntraLeadInterReal = fillfunc @@@ {{5;;, 6}, {5;;, 3}};
		Map[HBlochFull[vk, #]&, {HLeadIntraInterReal, HCSRIntraLeadInterReal}, {2}]),
		True, 0
	]
];*)

(*HBlochsForSpecFunc[vk_, ptscellsvas_, tfunc_, dup_] :=
Module[{HLeadIntraInterReal, HCSRIntraLeadInterReal, fillfunc, keys = Keys[ptscellsvas], len = Length[ptscellsvas]},
	(*fillfunc[indfs_, indi_] := HMatrixFromHoppings[{#, ptscellsvas[[indi]]}, tfunc, dup] & /@ KeyMap[# - keys[[indi]] &][ptscellsvas[[indfs]]];*)
	fillfunc[indfs_, indi_] := HMatrixFromHoppings[{#, ptscellsvas[[indi]]}, tfunc, dup] & /@ ptscellsvas[[indfs]];
	HLeadIntraInterReal = fillfunc @@@ {{2;;4, 3}, {2;;4, 1}};
	HCSRIntraLeadInterReal = fillfunc @@@ {{5;;, 6}, {5;;, 3}};
	Map[HBlochFull[vk, #]&, {HLeadIntraInterReal, HCSRIntraLeadInterReal}, {2}]
]*)

(*PhotonBlocks[{A0_, Avecn:(_Function|_Symbol), \[Omega]_}, mnup_Integer][ptf_, pti_] :=
Module[{vd = ptf - pti, zero = 1.*^-5, d, ele, dim = (2 mnup + 1){1, 1}, photondress, sparsezero, sparseid, sparsediag},
	(*A0 -> q A0/\[HBar], \[Omega] -> \[HBar] \[Omega]*)
	d = Norm[vd];
	ele[m_, n_] := 1/(2\[Pi]) NIntegrate[Exp[I (A0 vd . Avecn[\[CurlyPhi]] + (m - n) \[CurlyPhi])], {\[CurlyPhi], -\[Pi], \[Pi]}, Method -> "LocalAdaptive"] // Chop;
	(*ele[m_,n_]:=Module[{\[CurlyPhi]},
	\[CurlyPhi]=ArcTan@@Reverse[vd]+\[Pi];
	Exp[-\[ImaginaryI](m-n) \[CurlyPhi]]BesselJ[m-n,A0 d]//N
	];*)
	photondress := Array[ele, dim, -mnup] // Chop;
	(*sparseconst[i_] := ConstantArray[i, dim, SparseArray];(*\:4f4e\:7ea7\:9519\:8bef\:ff01\:ff01\:ff01*)*)
	sparsezero = ConstantArray[0, dim, SparseArray];
	sparseid = IdentityMatrix[dim, SparseArray];
	sparsediag := SparseArray[Band[{1, 1}]-> -\[Omega] Range[-mnup, mnup]];
	If[d > zero,
		{photondress, sparsezero},
		{sparseid, sparsediag}
	]
];*)

(*Options[PhotonBlocksTensor] = Options[FourierCoefficient];
PhotonBlocksTensor[functime: (_Function|_Symbol), \[Omega]_, mnup_Integer, opts:OptionsPattern[]][ptf_, pti_] :=
Module[{vd = ptf - pti, zero = 1.*^-5, d, coef, ele, dim = (2 mnup + 1){1, 1}, innerdof = Dimensions[functime[0.123]], photondress, sparseid, sparsezero, sparsediag, innerid},
	d = Norm[vd]; innerid = IdentityMatrix[innerdof, SparseArray];
	coef[l_Integer] := coef[l] = FourierCoefficient[functime[\[CurlyPhi]], \[CurlyPhi], l, opts] // FullSimplify;
	ele[m_, n_] := coef[n - m];
	photondress = Array[ele, dim, -mnup];
	sparsezero := ConstantArray[0, dim ~Join~ innerdof, SparseArray];
	sparseid := TensorProduct[IdentityMatrix[dim, SparseArray], innerid];
	sparsediag := TensorProduct[SparseArray[Band[{1, 1}]-> -\[Omega] Range[-mnup, mnup]], innerid];
	If[d > zero,
		{photondress, sparsezero},
		{sparseid, photondress + sparsediag}]
];

Options[NPhotonBlocksTensor] = Options[NIntegrate];
NPhotonBlocksTensor[functime: (_Function|_Symbol), \[Omega]_, mnup_Integer, opts:OptionsPattern[]][ptf_, pti_] :=
Module[{vd = ptf - pti, zero = 1.*^-5, d, coef, ele, dim = (2 mnup + 1){1, 1}, innerdof = Dimensions[functime[0.123]], photondress, sparseid, sparsezero, sparsediag, innerid},
	d = Norm[vd]; innerid = IdentityMatrix[innerdof, SparseArray];
	coef[l_Integer] := coef[l] = 1/(2\[Pi]) NIntegrate[functime[\[CurlyPhi]] Exp[-I l \[CurlyPhi]], {\[CurlyPhi], -\[Pi], \[Pi]}, opts, AccuracyGoal -> 10, Method -> "LocalAdaptive"] // Chop;
	ele[m_, n_] := coef[n - m];
	photondress = Array[ele, dim, -mnup] // Chop;
	sparsezero := ConstantArray[0, dim ~Join~ innerdof, SparseArray];
	sparseid := TensorProduct[IdentityMatrix[dim, SparseArray], innerid];
	sparsediag := TensorProduct[SparseArray[Band[{1, 1}]-> -\[Omega] Range[-mnup, mnup]], innerid];
	If[d > zero,
		{photondress, sparsezero},
		{sparseid, photondress + sparsediag}]
];

PhotonDressTensor[t_, photonblocks_] :=
Module[{n = Length[t], s, combine},
	s = IdentityMatrix[n, SparseArray];
	combine = TensorContract[TensorProduct[#2, #], {{3, 6}}] &;
	ArrayFlatten[MapThread[combine, {{t, s}, photonblocks}] // Total]
];*)

(*HFloquetEffectiveBlochMatrixFromExtended[\[Omega]_, mnup_Integer][hmatrixfromhopping_] :=
Module[{dim = 2mnup + 1, Hs, H0},
	(*Hs = Transpose[Partition[hmatrixfromhopping, dim{1, 1}], {3, 4, 1, 2}] // SparseArray;*)
	Hs = flqReprAlt[dim][hmatrixfromhopping];
	H0 = Hs[[mnup + 1, mnup + 1]];
	H0 + 1/\[Omega] Sum[comm[Hs[[1, i]], Hs[[i, 1]]]/(i - 1), {i, 2, dim}]
];*)
(*Options[HFloquetEffectiveBlochMatrixFromExtended] = {"ExpansionOrder" -> 1, "ReturnComponents" -> False};
HFloquetEffectiveBlochMatrixFromExtended[\[Omega]_, mnup_Integer, opts:OptionsPattern[]][hmatrixfromhopping_] :=
Module[{dim = 2mnup + 1, lmax = 2mnup, mnrange, Hs, H, H1, H2pure, H2mixed, components, order = OptionValue["ExpansionOrder"]},
	mnrange = DeleteCases[Range[-lmax, lmax], 0];
	Hs = flqReprAlt[dim][hmatrixfromhopping];
	H[0] = Hs[[mnup + 1, mnup + 1]];
	Do[{H[-(i-1)], H[i-1]} = {Hs[[1, i]], Hs[[i, 1]]}, {i, 2, dim}];
	H[l_Integer /; Abs[l] > lmax] := 0 H[0];
	
	H1 = Sum[comm[H[-l], H[l]]/l, {l, lmax}];
	H2pure := Sum[comm[H[-m], comm[H[0], H[m]]]/(2 m^2), {m, mnrange}];
	H2mixed := Sum[comm[H[-n], comm[H[n - m], H[m]]]/(3 m n), {m, mnrange}, {n, DeleteCases[mnrange, m]}];
	components = Which[
		order == 1, <|"Order0" -> H[0], "Order1" -> H1/\[Omega]|>,
		order == 2, <|"Order0" -> H[0], "Order1" -> H1/\[Omega], "Order2Pure" -> H2pure/\[Omega]^2, "Order2Mixed" -> H2mixed/\[Omega]^2|>
		];
	If[OptionValue["ReturnComponents"], components, Chop[Total[components]]]
];*)

(*HFloquetEffectiveHoppingMatricesFromExtended[\[Omega]_, mnup_Integer][h0isvas_Association] :=
Module[{dim = 2mnup + 1, hsvas, hefffunc, H0, h0s0i, h0isvasrepralt, h0isvasalleffective, zero = 1.*^-5},
	hefffunc[{va1_ -> h1_, va2_ -> h2_}] := (va1 + va2) -> 1/\[Omega] Sum[comm[h1[[1, i]], h2[[i, 1]]]/(i - 1), {i, 2, dim}];
	(*h0isvasrepralt = SparseArray[Transpose[Partition[#, dim{1, 1}], {3, 4, 1, 2}]] & /@ h0isvas;*)
	h0isvasrepralt = flqReprAlt[dim] /@ h0isvas;
	h0s0i = Normal[h0isvasrepralt[[;;, mnup + 1, mnup + 1]]];
	hsvas = Normal[h0isvasrepralt];
	h0isvasalleffective = Join[h0s0i, hefffunc /@ Tuples[hsvas, 2]];
	DeleteCases[s_ /; s["Density"] < zero][GroupBy[h0isvasalleffective, Keys -> Values, Chop @* Total]]
];*)
(*Options[HFloquetEffectiveHoppingMatricesFromExtended] = {"ExpansionOrder" -> 1, "ReturnComponents" -> False};
HFloquetEffectiveHoppingMatricesFromExtended[\[Omega]_, mnup_Integer, opts : OptionsPattern[]][h0isvas_Association] :=
Module[{dim = 2 mnup + 1, lmax = 2 mnup, mnrange, h0isvasrepralt, hsvas, h0s0i, h, resultant, pairgroups, triplegroups, hefffunc1, hefffunc2pure,
		hefffunc2mixed, H0, H1, H2pure, H2mixed, components, order = OptionValue["ExpansionOrder"], zero = 1.*^-5, clean},
    h0isvasrepralt = flqReprAlt[dim] /@ h0isvas;
    hsvas = Normal[h0isvasrepralt];
    h0s0i = Normal[h0isvasrepralt[[;;, mnup + 1, mnup + 1]]];
    h[hmat_, l_Integer] := Which[l == 0, hmat[[mnup + 1, mnup + 1]],
            1 <= l <= lmax, hmat[[l + 1, 1]],
	        -lmax <= l <= -1, hmat[[1, 1 - l]],
			True, 0 hmat[[mnup + 1, mnup + 1]]];
    resultant[tuple_List] := Total[First /@ tuple];
    clean[assoc_(*Association*)] := DeleteCases[s_ /; s["Density"] < zero][assoc];
    mnrange = DeleteCases[Range[-lmax, lmax], 0];
    hefffunc1[{_ -> h1_, _ -> h2_}] := 1/\[Omega] Sum[comm[h[h1, -l], h[h2, l]]/l, {l, 1, lmax}];
    hefffunc2pure[{_ -> h1_, _ -> h2_, _ -> h3_}] := 1/\[Omega]^2 Sum[comm[h[h1, -m], comm[h[h2, 0], h[h3, m]]]/(2 m^2), {m, mnrange}];
    hefffunc2mixed[{_ -> h1_, _ -> h2_, _ -> h3_}] := 1/\[Omega]^2 Sum[comm[h[h1, -n], comm[h[h2, n - m], h[h3, m]]]/(3 m n), {m, mnrange}, {n, DeleteCases[mnrange, m]}];
    H0 = clean[h0s0i];
    pairgroups = GroupBy[Tuples[hsvas, 2], resultant];
    H1 = clean[Chop[Total[hefffunc1 /@ #]] & /@ pairgroups];
    If[order == 2, triplegroups = GroupBy[Tuples[hsvas, 3], resultant];
        H2pure = clean[Chop[Total[hefffunc2pure /@ #]] & /@ triplegroups];
        H2mixed = clean[Chop[Total[hefffunc2mixed /@ #]] & /@ triplegroups];
    ];
    components = Which[
        order == 1, <|"Order0" -> H0, "Order1" -> H1|>,
        order == 2, <|"Order0" -> H0, "Order1" -> H1, "Order2Pure" -> H2pure, "Order2Mixed" -> H2mixed|>
    ];
    If[OptionValue["ReturnComponents"], components,
        clean[Merge[Values[components], Chop @* Total]]]
];*)

(*HCSRBlocksAndersonDisordered[W_, nensemble_Integer:1, innerdofmat_][hdods:{{__}, {__}}] :=
Module[{hds, hods, hdsdisordered, hdsensemble, innerdof = Length[innerdofmat], disorderblock},
	{hds, hods} = hdods;
	disorderblock = DiagonalMatrix[RandomReal[{-1, 1} W/2, Length[#]/innerdof], TargetStructure -> "Sparse"] &;
	hdsdisordered := # + KroneckerProduct[disorderblock[#], innerdofmat] & /@ hds;
	Table[{hdsdisordered, hods}, nensemble]
];*)

(*\[Lambda]SKOnsiteBlock["sp"][\[Lambda]_] := ConstantArray[0, {1, 3}2, SparseArray];
\[Lambda]SKOnsiteBlock["sd"][\[Lambda]_] := ConstantArray[0, {1, 5}2, SparseArray];
\[Lambda]SKOnsiteBlock["pd"][\[Lambda]_] := ConstantArray[0, {3, 5}2, SparseArray];*)




(* ::Title::Closed:: *)
(*EigenSpect`**)


(*FirstBrillouinZone::usage = "Shows the first Brillouin zone with reciprocal lattice vectors.";*)
(*PlaquetteChern::usage = "Calculates the Chern number contributed by a small region, usually dubbed as a plaquette; summed over the first Brillouin zone cover, the full Chern number is obtained.";
PlaquetteRegionPartition::usage = "Partitions a region, e.g. the first Brillouin zone, with a triangularization cover.";*)

(*BandData[hbloch_, ks_, map_: Map, s:OptionsPattern[Eigenvalues]] :=
(Sort @ Eigenvalues[hbloch[#], s, Method -> "Direct"] & ~map~ ks)\[Transpose];*)
(*funcpattern = (_Function | _Symbol | _[__]);*)
(*funcpattern = (_Function | _Symbol);*)

(*Options[ParallelBandDataWithWeight] = Join[Options[Eigensystem], {"StateFunction" -> (Total[Abs[#]^4] &)}];
ParallelBandDataWithWeight[h_, kgrid_, n_Integer, ps:OptionsPattern[]] :=
Module[{ps2},
	ps2 = Sequence @@ FilterRules[{ps}, Options[Eigensystem]];
	MapAt[OptionValue["StateFunction"], Eigensystem[h[#], n, ps2, Method -> {"Arnoldi", "MaxIterations" -> \[Infinity]}]\[Transpose] // Sort, {All, 2}] & ~ParallelMap~ kgrid
];
(*parallelBandDataWithState[h_,kgrid_,n_Integer,infoextracfunc_:(Total[Abs[#]^4]^(1/4)&),ps:OptionsPattern[Eigensystem]]:=MapAt[infoextracfunc,Eigensystem[h[#],n,ps,Method->{"Arnoldi","MaxIterations"->\[Infinity]}(*Method\[Rule]"Banded"*)]\[Transpose]//Sort,{All,2}]&~ParallelMap~kgrid;*)

Options[BandDataWithWeight] = Join[Options[Eigensystem], {"StateFunction" -> (Total[Abs[#]^4] &)}];
BandDataWithWeight[h_, kgrid_, n_Integer, ps:OptionsPattern[]]:=
Module[{ps2},
	ps2 = Sequence @@ FilterRules[{ps}, Options[Eigensystem]];
	MapAt[OptionValue["StateFunction"], Eigensystem[h[#], n, ps2, Method -> {"Arnoldi", "MaxIterations" -> \[Infinity]}]\[Transpose] // Sort, {All, 2}] & ~Map~ kgrid
];*)

(*ParallelBandDataWithWeight[h_, kgrid_, n_Integer, ps:OptionsPattern[]] :=
Module[{ps1val, ps2, eigensyst},
	ps1val = OptionValue["StateFunction"];
	ps2 = TBMethod`DataVisualization`Private`optionsselect[ps, Eigensystem];
	eigensyst := Sort[Eigensystem[h[#], n, ps2, Method -> {"Arnoldi", "MaxIterations" -> \[Infinity]}]\[Transpose]] & ~ParallelMap~ kgrid;
	Which[
		MatchQ[ps1val, _Function], MapAt[ps1val, {All, All, 2}][eigensyst],
		MatchQ[ps1val, {__Function}], Transpose[Map[Thread, MapAt[Through @* ps1val, {All, All, 2}][eigensyst], {2}], {2, 3, 1}],
		True, Message[ParallelBandDataWithWeight::weightfunc, ps1val]
	]
];*)

(*BandDataWithWeight[h_, kgrid_, n_Integer, ps:OptionsPattern[]]:=
Module[{ps1val, ps2, eigensyst},
	ps1val = OptionValue["StateFunction"];
	ps2 = TBMethod`DataVisualization`Private`optionsselect[ps, Eigensystem];
	eigensyst := Sort[Eigensystem[h[#], n, ps2, Method -> {"Arnoldi", "MaxIterations" -> \[Infinity]}]\[Transpose]] & ~Map~ kgrid;
	Which[
		MatchQ[ps1val, _Function], MapAt[ps1val, {All, All, 2}][eigensyst],
		MatchQ[ps1val, {__Function}], Transpose[Map[Thread, MapAt[Through @* ps1val, {All, All, 2}][eigensyst], {2}], {2, 3, 1}],
		True, Message[BandDataWithWeight::weightfunc, ps1val]
	]
];*)

(*FirstBrillouinZone[vbs_ /; Dimensions[vbs] == {2, 2} || Dimensions[vbs] == {3, 3} , n_:2, opts:OptionsPattern[Show]] :=
Module[{intcoeffs, fbz, reciprocalvectors, len = Length[vbs], \[CapitalGamma], colors},
	intcoeffs = Tuples[Range[-n, n], len];
	\[CapitalGamma] = ConstantArray[0, len];
	colors = Take[{Red, Green, Blue}, len];
	fbz = MinimalBy[Norm @* RegionCentroid][MeshPrimitives[VoronoiMesh[intcoeffs . vbs], len]] // First;
	reciprocalvectors = MapThread[{#, Arrow[{\[CapitalGamma], #2}]} &, {colors, vbs}];
	Show[{HighlightMesh[fbz, {Style[1, Black, Thick], Style[2, Gray, Opacity[0.1]]}],
	If[len == 2, Graphics, Graphics3D][{Thick, reciprocalvectors}]}, opts]
];*)

(*FirstBrillouinZoneRegion[vbs_ /; Dimensions[vbs] == {2, 2} || Dimensions[vbs] == {3, 3} , n_:2] :=
Module[{intcoeffs, reciprocalvectors, len = Length[vbs]},
	intcoeffs = Tuples[Range[-n, n], len];
	MinimalBy[Norm @* RegionCentroid][MeshPrimitives[VoronoiMesh[intcoeffs . vbs], len]] // First
];(*suffers RAM limit suddenly*)*)

(*LabelPathSamplings[pathsamplings_, labels:{__String}] :=
Module[{func, lbllen = Length[labels], numberlen, ptsall, numbers, numbersfinal},
	func = MapAt[lis |-> Callout[lis, #2[[1]], Automatic, Automatic, Appearance -> "CurvedLeader"], #2[[2]]][#] &;
	{ptsall, numbers} = pathsamplings; numberlen = Length[numbers];
	numbersfinal = If[lbllen == numberlen-1, Most @ numbers, MapAt[#-1 &, -1] @ numbers];
	Fold[func, ptsall, {labels, numbersfinal}\[Transpose]]
];*)

(*berryCurvatureCore[
	{Hveck_?MatrixQ, Vveck : {__?MatrixQ}}, Q_, nF_?(# \[Element] PositiveIntegers &),
	dirs : {\[Alpha]_Integer, \[Beta]_Integer}, opts:OptionsPattern[Eigensystem]
] /; 1 <= \[Alpha] <= Length[Vveck] && 1 <= \[Beta] <= Length[Vveck] && (*\[Alpha] != \[Beta] && *)nF < Length[Hveck] :=
Module[{q, v\[Alpha], v\[Beta], j\[Alpha], eigensyst, \[CapitalOmega]nl},
	q =If[MatrixQ[Q], Q, Q IdentityMatrix[Length[Hveck], SparseArray]];
	{v\[Alpha], v\[Beta]} = Vveck[[dirs]]; j\[Alpha] = (q . v\[Alpha] + v\[Alpha] . q)/2;
	eigensyst = Sort[Eigensystem[Hveck, opts, Method -> "Direct"]\[Transpose]];
	\[CapitalOmega]nl = (#2[[2]]\[Conjugate] . j\[Alpha] . #[[2]]	#[[2]]\[Conjugate] . v\[Beta] . #2[[2]])/(#2[[1]] - #[[1]])^2 &;
	2 Im @ Total[\[CapitalOmega]nl @@@ Tuples[TakeDrop[eigensyst, nF]]]
];
BerryCurvature[
	{H_, V_}, Q_, nF_?(# \[Element] PositiveIntegers &), 
	dirs : {_Integer, _Integer} : {1, 2},
	opts:OptionsPattern[Eigensystem]
][veck_List] := berryCurvatureCore[{H[veck], V[veck]}, Q, nF, dirs, opts];
BerryCurvature[
	{H_, V_}, nF_?(# \[Element] PositiveIntegers &),
	dirs : {_Integer, _Integer} : {1, 2},
	opts : OptionsPattern[Eigensystem]
][veck_List] := berryCurvatureCore[{H[veck], V[veck]}, 1, nF, dirs, opts];
BerryCurvature[
	HV_, Q_, nF_?(# \[Element] PositiveIntegers &), 
	dirs : {_Integer, _Integer} : {1, 2},
	opts:OptionsPattern[Eigensystem]
][veck_List] := berryCurvatureCore[HV[veck], Q, nF, dirs, opts];
BerryCurvature[
	HV_, nF_?(# \[Element] PositiveIntegers &),
	dirs : {_Integer, _Integer} : {1, 2},
	opts : OptionsPattern[Eigensystem]
][veck_List] := berryCurvatureCore[HV[veck], 1, nF, dirs, opts];*)

(*quantumGeometricCore[
	{Hveck_?MatrixQ, Vveck : {__?MatrixQ}}, Q_,
	nF_?(# \[Element] PositiveIntegers &),
	dirs : {\[Alpha]_Integer, \[Beta]_Integer},
	opts : OptionsPattern[Eigensystem]
] /; 1 <= \[Alpha] <= Length[Vveck] && 1 <= \[Beta] <= Length[Vveck] && (*\[Alpha] != \[Beta] && *)nF < Length[Hveck] && (NumericQ[Q] || (MatrixQ[Q] && Dimensions[Q] === Dimensions[Hveck])) :=
Module[{q, v\[Alpha], v\[Beta], j\[Alpha], eigensyst, \[Chi]nl},
	q = If[MatrixQ[Q], Q, Q IdentityMatrix[Length[Hveck], SparseArray]];
	{v\[Alpha], v\[Beta]} = Vveck[[dirs]];
	j\[Alpha] = (q . v\[Alpha] + v\[Alpha] . q)/2;
	eigensyst = Sort[Eigensystem[Hveck, opts, Method -> "Direct"]\[Transpose]];
	\[Chi]nl = (#2[[2]]\[Conjugate] . j\[Alpha] . #[[2]] #[[2]]\[Conjugate] . v\[Beta] . #2[[2]]) / (#2[[1]] - #[[1]])^2 &;
	Total[\[Chi]nl @@@ Tuples[TakeDrop[eigensyst, nF]]]
];*)

(*WannerChargeCenter[] :=.*)

(*PlaquetteChern[vks:{{__?NumericQ}..}, heff_, nF_?(# \[Element] PositiveIntegers &), opts:OptionsPattern[Eigensystem]] :=
Module[{occupiedstates, matFs, stateloop, func},
	occupiedstates = Take[Sort[Eigensystem[heff[#], opts, Method -> "Direct"]\[Transpose]], nF][[;;, 2]] & /@ vks;
	func = {vecs1, vecs2} |-> Outer[#\[Conjugate] . #2 &, vecs1, vecs2, 1];
	stateloop = Partition[occupiedstates, 2, 1, {1, 1}];
	matFs = Dot @@ func @@@ stateloop;
	(*-1/(2\[Pi]) Arg[Eigenvalues[matFs]] // Sort*)
	-(1/(2\[Pi])) Arg @ Det[matFs]
];

PlaquetteRegionPartition[region_, opts:OptionsPattern[TriangulateMesh]] :=
Module[{regiondiscrized, meshcoordinates, plaquettevertexindex},
	regiondiscrized = TriangulateMesh[region, opts];
	meshcoordinates = MeshCoordinates[regiondiscrized];
	Echo[MeshRegion[regiondiscrized, PlotTheme -> "Lines", PlotLabel -> StringTemplate["Vertex number: ``"][Length[meshcoordinates]]]];
	plaquettevertexindex = MeshCells[regiondiscrized, 2][[;;, 1]];
	Extract[meshcoordinates, {#}\[Transpose]] & /@ plaquettevertexindex
];*)

(*PlaquetteRegionPartitionComplex[region_, opts:OptionsPattern[TriangulateMesh]] :=
Module[{regiondiscrized, meshcoordinates, plaquettevertexindex},
	regiondiscrized = TriangulateMesh[region, opts, 
		MaxCellMeasure -> {"Area" -> (Area[region].01)}, Method -> "ConstrainedQuality"];
	meshcoordinates = MeshCoordinates[regiondiscrized];
	plaquettevertexindex = List @@@ MeshCells[regiondiscrized, 2];
	Echo[
		MeshRegion[regiondiscrized, PlotTheme -> "Lines", 
		PlotLabel -> StringTemplate["Vertex #: ``, Plaquette #: ``."][Length[meshcoordinates], Length[plaquettevertexindex]]]
	];
	{meshcoordinates, plaquettevertexindex}
];*)

(*plaquettePhase[occupiedstates_] :=
Module[{stateloop, matD, func},
	stateloop = Partition[occupiedstates, 2, 1, {1, 1}];
	func = {vecs1, vecs2} |-> Outer[#\[Conjugate] . #2 &, vecs1, vecs2, 1];
	matD = Dot @@ func @@@ stateloop;
	Arg @ Det[matD]
];*)

(*plaquettePhase[occupiedstates_] :=
Module[{stateloop, matD, func},
	stateloop = {#, RotateLeft[#]} &[occupiedstates];
	func = {vecs1, vecs2} |-> Outer[#\[Conjugate] . #2 &, vecs1, vecs2, 1];
	matD = Dot @@ MapThread[func][stateloop];
	Arg @ Det[matD]
];*)




(* ::Title::Closed:: *)
(*LGFF`**)


(* The transfer matrix *)
(*Tmat[matlist_] := Module[{positions(*, n = Length[matlist]*)},
	positions[m_] := Normal @ Partition[SparseArray[{{1} -> 1, {i_?OddQ} :> (i + 1)/2, {-1} -> 1}, 2 m, 2], 2];
	SparseArray[ Total[ Dot @@ Extract[matlist, positions[#] ] & /@ Range[Length[matlist]] ] ]
];*)

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

(*SurfaceGreen[e_, {h0_, h1_}, mode:(1|2|3):3] :=
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
];*)

(*Sigma[epsilon_, {h0_, h1_, H01_}, mode:(1|2|3):1] :=
Module[{gsurface},
	gsurface = SurfaceGreen[epsilon, {h0, h1}, mode];
	SparseArray[H01 . gsurface . H01\[ConjugateTranspose]]
];*)

(*LocalDOSRealSpace[Gs_, sigmas_, layeredpts_Association, innerdof_:1] :=
Module[{gamma, innerdofldos, ldos},
	gamma = I (# - ConjugateTranspose[#]) & @ Total[sigmas];
	(*gamma = Im @ Total[sigmas];*)
	innerdofldos = (*-*)1/\[Pi] Diagonal[ # . gamma . ConjugateTranspose[#] ] & /@ Gs;
	(*innerdofldos = -1/\[Pi] Diagonal[ Im[#] ] & /@ Gs;*) (*WRONG!*)
	ldos = BlockMap[Total, #, innerdof] & /@ innerdofldos;
	MapThread[Append, Join @@@ {Values[layeredpts], Reverse[ldos]}]
];*)

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
(*LocalDOSReciprocalSpace[\[Epsilon]_, {HLeadBloch_?MatrixQ, HLead12_?MatrixQ}, mode:(1|2|3):3] := LocalDOSReciprocalSpace[\[Epsilon], {{HLeadBloch, HLead12}, {HLeadBloch, HLead12}}, mode];
LocalDOSReciprocalSpace[\[Epsilon]_, {HLeadBloch_?MatrixQ, HLead12_?MatrixQ}, HCSRBloch_?MatrixQ, mode:(1|2|3):3] := LocalDOSReciprocalSpace[\[Epsilon], {{HLeadBloch, HLead12}, {HCSRBloch, HLead12}}, mode];
LocalDOSReciprocalSpace[\[Epsilon]_, {{HLeadBloch_?MatrixQ, HLead12_?MatrixQ}, {HCSRBloch_?MatrixQ, HCSRLead1_?MatrixQ}}, mode:(1|2|3):3] :=
Module[{zero = 1.*^-4, \[CapitalSigma]},
	\[CapitalSigma] = Sigma[Complex[\[Epsilon], zero], {HLeadBloch, HLead12, HCSRLead1}, mode];
	-Im @ Tr @ CentralGreen[Complex[\[Epsilon], zero], HCSRBloch, {\[CapitalSigma]}]
];*)

(*currentTensorBlocks[Gs_, sigmas_, blockHs: {ds_, os_}, innerdof_:1] :=
Module[{gamma, blockGns, Gsre = Reverse[Gs], jblock0, jblock0innersummed},
	gamma = I (# - #\[ConjugateTranspose]) & [Total[sigmas]];
	blockGns = {# . gamma . #\[ConjugateTranspose] & /@ Gsre, # . gamma . #2\[ConjugateTranspose] & @@@ Partition[Gsre, 2, 1]};
	(*jblock0 = -Im[MapAt[ConjugateTranspose, {2, All}][blockHs] blockGns];*)
	(*jblock0 = Im[MapAt[Transpose, {2, All}][blockHs] blockGns];*)
	jblock0 = Im[Map[Transpose, blockHs, {2}] blockGns];(*!!!*)
	jblock0innersummed = Table[BlockMap[Total[#, 2] &, #, {1, 1}innerdof] & /@ x, {x, jblock0}];
	(*How to sum the internal degree of freedom?*)
	Append[jblock0innersummed, -Transpose /@ jblock0innersummed[[2]]]
];(*bond current in layered block form*)*)
(*currentTensorBlocks[innerproj_][Gs_, sigmas_, blockHs: {ds_, os_}] :=
Module[{gamma, blockGns, Gsre = Reverse[Gs], jblock0, jblock0innersummed, innerdof = Dimensions[innerproj]},
	gamma = I (# - #\[ConjugateTranspose]) & [Total[sigmas]];
	blockGns = {# . gamma . #\[ConjugateTranspose] & /@ Gsre, # . gamma . #2\[ConjugateTranspose] & @@@ Partition[Gsre, 2, 1]};
	(*jblock0 = -Im[MapAt[ConjugateTranspose, {2, All}][blockHs] blockGns];*)
	(*jblock0 = Im[MapAt[Transpose, {2, All}][blockHs] blockGns];*)
	jblock0 = Im[Map[Transpose, blockHs, {2}] blockGns];(*!!!*)
	jblock0innersummed = Table[BlockMap[Total[innerproj . # . innerproj, 2] &, #, {1, 1}innerdof] & /@ x, {x, jblock0}];
	(*How to sum the internal degree of freedom?*)
	Append[jblock0innersummed, -Transpose /@ jblock0innersummed[[2]]]
];(*bond current in layered block form*)*)
(*currentTensorBlocks[proj_][Gs_, sigmas_, blockHs: {ds_, os_}] :=
Module[{gamma, blockGns, Gsre = Reverse[Gs], jblock0, jblock0innersummed, innerdof = Dimensions[proj], blockGnspartial, blockHspartial},
	gamma = I (# - #\[ConjugateTranspose]) & [Total[sigmas]];
	blockGns = {# . gamma . #\[ConjugateTranspose] & /@ Gsre, # . gamma . #2\[ConjugateTranspose] & @@@ Partition[Gsre, 2, 1]};
	blockGnspartial = Table[ArrayFlatten[BlockMap[proj . # &, #, innerdof]] & /@ x, {x, blockGns}];
	blockHspartial = Map[Transpose, blockHs, {2}];
	jblock0 = Im[blockHspartial blockGnspartial];
	jblock0innersummed = Table[BlockMap[Total[#, 2] &, #, innerdof] & /@ x, {x, jblock0}];
	(*How to sum the internal degree of freedom?*)
	Append[jblock0innersummed, -Transpose /@ jblock0innersummed[[2]]]
];(*bond current in layered block form*)*)

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

(*LocalCDV[Gs_, sigmas_, blockHs: {ds_, os_}, layeredpts_Association, innerdof_:1] :=
Module[{ptscsr = Join @@ layeredpts, innersummed, js, currenttensorblocks, jtensorfull},
	currenttensorblocks = currentTensorBlocks[Gs, sigmas, blockHs, innerdof];
	jtensorfull = jtensorFromBlocks[currenttensorblocks];
	js = jvecfield[ptscsr, jtensorfull];
	(*KeyValueMap[List] @ (Total /@ GroupBy[js, First -> Last])*)
	KeyValueMap[List] @* Merge[Total] @ js
];*)

(*Options[transmissionsfunc] = {"SigmaMode" -> 3, "LeadGaugeTransforms" -> None};
transmissionsfunc::glen = "\"LeadGaugeTransforms\" contains `1` transformations, but `2` leads are present.";
transmissionsfunc[\[Epsilon]_, hcsrdod_, leadshs_, opts : OptionsPattern[]] :=
Module[{\[CapitalSigma]s, \[CapitalSigma]s0, blockG, ter = Length[leadshs], gUs = OptionValue["LeadGaugeTransforms"]},
	If[gUs =!= None && (!ListQ[gUs] || Length[gUs] =!= ter),
		Message[transmissionsfunc::glen, If[ListQ[gUs], Length[gUs], "non-list"], ter];
        Return[$Failed]
    ];
	
	\[CapitalSigma]s0 = Sigma[\[Epsilon], #, OptionValue["SigmaMode"]] & /@ leadshs;
	\[CapitalSigma]s = If[gUs === None, \[CapitalSigma]s0, MapThread[# . #2 . (#)\[ConjugateTranspose] &, {gUs, \[CapitalSigma]s0}]];
    blockG = CentralBlockGreens[\[Epsilon], hcsrdod, \[CapitalSigma]s, "T"];
    Table[If[p == q, 0., Transmission[blockG, \[CapitalSigma]s[[{p, q}]]]], {p, ter}, {q, ter - 1}]
];*)

(*Options[transmissionsfunc] = {"SigmaMode" -> 3};
(*truncated transmission matrix for a grounded last terminal*)
transmissionsfunc[\[Epsilon]_, hcsrdod_, leadshs_, gUs_, opts:OptionsPattern[]]:=
Module[{\[CapitalSigma]s, blockG, \[CapitalSigma]s0, ter = Length[leadshs]},
	\[CapitalSigma]s0 = Sigma[\[Epsilon], #, OptionValue["SigmaMode"]] & /@ leadshs;
	\[CapitalSigma]s = MapThread[# . #2 . (#\[ConjugateTranspose]) &, {gUs, \[CapitalSigma]s0}];
	blockG = CentralBlockGreens[\[Epsilon], hcsrdod, \[CapitalSigma]s, "T"];
	(*Table[If[p == q || q == ter, 0., Transmission[blockG, \[CapitalSigma]s[[{p, q}]]]], {p, ter}, {q, ter}]*)
	Table[If[p == q, 0., Transmission[blockG, \[CapitalSigma]s[[{p, q}]]]], {p, ter}, {q, ter - 1}]
];

Options[HallAndLongitudinalResistances] = Options[transmissionsfunc];
HallAndLongitudinalResistances[\[Epsilon]_, hcsrdod_, leadshs_, gUs_, opts:OptionsPattern[]] :=
Module[{\[ScriptCapitalT], \[ScriptCapitalT]func, Rfunc, cnup = 1.*^7, ter = Length[leadshs], transmissions},
	(*\[ScriptCapitalT]func = DiagonalMatrix[Total[#]] - # &;*)
	\[ScriptCapitalT]func = DiagonalMatrix[Total[#]] - Most[#] &;
	Rfunc = # - {#2, #3} & @@ Rest[LinearSolve[#][UnitVector[ter - 1, 1]]] &;
	transmissions = transmissionsfunc[\[Epsilon], hcsrdod, leadshs, gUs, opts];
	(*\[ScriptCapitalT] = Drop[\[ScriptCapitalT]func[transmissions], -1, -1];*)
	\[ScriptCapitalT] = \[ScriptCapitalT]func[transmissions];
	If[LUDecomposition[\[ScriptCapitalT]][[4]] > cnup, {"NaN", "NaN"}, (*condition number from LU*)
		Rfunc[\[ScriptCapitalT]]
	]
];*)

(*Options[HallAndLongitudinalResistances] = {"SigmaMode" -> 3};
HallAndLongitudinalResistances[\[Epsilon]_, hcsrdod_, leadshs_, gUs_, opts:OptionsPattern[]] :=
Module[{\[ScriptCapitalT], \[ScriptCapitalT]func, Rfunc, cnup = 1.*^7, ter = Length[leadshs], transmissions},
	(*\[ScriptCapitalT]func = # - DiagonalMatrix[Total[#, {2}]] &;*)
	\[ScriptCapitalT]func = DiagonalMatrix[Total[#]] - # &;
	Rfunc = # - {#2, #3} & @@ Rest[LinearSolve[#][UnitVector[ter - 1, 1]]] &;
	transmissions = Module[{\[CapitalSigma]s, blockG, \[CapitalSigma]s0},
		\[CapitalSigma]s0 = Sigma[\[Epsilon], #, OptionValue["SigmaMode"]] & /@ leadshs;
		\[CapitalSigma]s = MapThread[# . #2 . (#\[ConjugateTranspose]) &, {gUs, \[CapitalSigma]s0}];
		blockG = CentralBlockGreens[\[Epsilon], hcsrdod, \[CapitalSigma]s, "T"];
		(*Table[If[p == q || p == ter, 0., Transmission[blockG, \[CapitalSigma]s[[{p, q}]]]], {p, ter}, {q, ter}]*)
		Table[If[p == q || q == ter, 0., Transmission[blockG, \[CapitalSigma]s[[{p, q}]]]], {p, ter}, {q, ter}]
	];
	(*\[ScriptCapitalT] = Drop[\[ScriptCapitalT]func[-transmissions], -1, -1];*)
	\[ScriptCapitalT] = Drop[\[ScriptCapitalT]func[transmissions], -1, -1];
	(*If[LinearAlgebra`Private`MatrixConditionNumber[\[ScriptCapitalT]] > cnup, {"NaN", "NaN"},
		Rfunc[\[ScriptCapitalT]]
	]*)
	If[LUDecomposition[\[ScriptCapitalT]][[(*3*)4]] > cnup, {"NaN", "NaN"}, (*condition number from LU*)
		Rfunc[\[ScriptCapitalT]]
	]
];*)

(*Options[HallAndLongitudinalConductances] = Options[HallAndLongitudinalResistances];
HallAndLongitudinalConductances[\[Epsilon]_, hcsrdod_, leadshs_, gUs_, opts:OptionsPattern[]] :=
Module[{RH, RL},
	{RH, RL} = HallAndLongitudinalResistances[\[Epsilon], hcsrdod, leadshs, gUs, opts];
	{-RH, RL}/(RH^2 + RL^2)
];*)

(*HallAndLongitudinalResistances[\[Epsilon]_, hcsrdod_, leadshs_] :=
Module[{\[ScriptCapitalT], \[ScriptCapitalT]func, Rfunc, cnup = 1.*^7, ter = Length[leadshs], transmissions},
	\[ScriptCapitalT]func = # - DiagonalMatrix[Total[#, {2}]] &;
	Rfunc = # - {#2, #3} & @@ Rest[LinearSolve[#][UnitVector[ter - 1, 1]]] &;
	transmissions = Module[{\[CapitalSigma]s, blockG},
		\[CapitalSigma]s = Sigma[\[Epsilon], #, 3] & /@ leadshs;
		blockG = CentralBlockGreens[\[Epsilon], hcsrdod, \[CapitalSigma]s, "T"];
		Table[If[p == q || p == ter, 0., Transmission[blockG, \[CapitalSigma]s[[{p, q}]]]], {p, ter}, {q, ter}]
	];
	\[ScriptCapitalT] = Drop[\[ScriptCapitalT]func[-transmissions], -1, -1];
	(*If[LinearAlgebra`Private`MatrixConditionNumber[\[ScriptCapitalT]] > cnup, {"NaN", "NaN"},
		Rfunc[\[ScriptCapitalT]]
	]*)
	If[LUDecomposition[\[ScriptCapitalT]][[3]] > cnup, {"NaN", "NaN"}, (*condition number from LU*)
		Rfunc[\[ScriptCapitalT]]
	]
];

HallAndLongitudinalConductances[\[Epsilon]_, hcsrdod_, leadshs_] :=
Module[{RH, RL},
	{RH, RL} = HallAndLongitudinalResistances[\[Epsilon], hcsrdod, leadshs];
	{RH, RL}/(RH^2 + RL^2)
];*)

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


(* ::Title::Closed:: *)
(*DataVisalization`**)


(*RealSpaceLocalDOSPlot[evalandevec_List, ptsdisk:{{_, _, _}..}|{{_, _}..}, region_?RegionQ, innerdof_Integer:1, ops:OptionsPattern[Graphics]] :=
Module[{op, largestcomps, \[Eta] = 1.*^-4, len = Length[ptsdisk], ratio = 2},
	op = If[MatchQ[ptsdisk, {{_, _, _}..}], KeyValueMap[Append] @* (data |-> GroupBy[data, (#[[;;2]] &) -> Last, Total]), Identity];
	largestcomps = Select[Last[#] > \[Eta] &] @ Join[ptsdisk, {BlockMap[Total, Abs[evalandevec[[2]]]^2, innerdof]}\[Transpose], 2];
	Graphics[
	{{Opacity[.1], Green, region},
	 {Opacity[.4], Red, Disk[{#, #2}, Sqrt[len]ratio #3] & @@@ op[largestcomps]}},
	ops, PlotLabel -> StringTemplate["\!\(\*SubscriptBox[\(E\), \(\[VeryThinSpace]\)]\) = ``"][evalandevec[[1]]]]
] /; (Length[Partition[evalandevec[[2]], innerdof]] == Length[ptsdisk]);*)
(*RealSpaceLocalDOSPlot[evalandevec_List, ptsdisk:{{_, _, _}..}|{{_, _}..}, region_?RegionQ, innerdof_Integer:1, ratio_:2, ops:OptionsPattern[Graphics]] :=
Module[{op, largestcomps, \[Eta] = 1.*^-4, len = Length[ptsdisk], (*ratio = 2,*) evaldisp},
	op = If[MatchQ[ptsdisk, {{_, _, _}..}], KeyValueMap[Append] @* (data |-> GroupBy[data, (#[[;;2]] &) -> Last, Total]), Identity];
	largestcomps = Select[Last[#] > \[Eta] &] @ Join[ptsdisk, {BlockMap[Total, Abs[evalandevec[[2]]]^2, innerdof]}\[Transpose], 2];
	evaldisp = ToString[ScientificForm[Re @ evalandevec[[1]], 4], StandardForm];
	Graphics[{
		{FaceForm[{Opacity[.2], Green}], EdgeForm[Black], region},
		{Opacity[.4], Red, Disk[{#, #2}, Sqrt[len]ratio #3] & @@@ op[largestcomps]},
		{Text[StringTemplate["\!\(\*SubscriptBox[\(E\), \(\[VeryThinSpace]\)]\) = ``"][evaldisp]]}
		},
		ops
	]
] /; (Length[Partition[evalandevec[[2]], innerdof]] == Length[ptsdisk]);*)

(*LocalDOSTidy[data_, quantile_] :=
Module[{maxquant, clipped, min = Min[data]},
	maxquant = Quantile[data // Flatten, quantile];
	clipped = Clip[data, {min, maxquant}];
	GaussianFilter[clipped, 2]
];*)

(*Options[BandPlotWithWeight] = Join[Options[Graphics], Options[BarLegend], {Joined -> True, ColorFunction -> (Hue[2(1 - #)/3] &)}];
(*bandPlotWithWeight[banddatawithstate_,cfunc_,cname_String,joined_:(True|False),ps:OptionsPattern[Graphics]]:=*)
BandPlotWithWeight[banddatawithweight_,
				   hisymmptname : {(_String|OverBar[_String])..} : {""},
				   ptsnumbers : {_?NumericQ..} : {1},
				   yticks :{{_, _}..} : Automatic,
				   ps:OptionsPattern[]] :=
Module[{kbdat, colors, m, n, lines, bfig, legend, fontfamily = (*"Helvetica"*)(*"Times New Roman"*)"Arial", style,
		style2, bdat, cdat, cfunc = OptionValue[ColorFunction], dticks, frameticks, ps1, ps2},
	{bdat, cdat} = Transpose[banddatawithweight, {3, 2, 1}]; {m, n} = Dimensions[bdat];
	dticks = {ptsnumbers, hisymmptname}\[Transpose]; frameticks = {{(*Automatic*)yticks, None}, {dticks, None}};
	style = {FontSize -> 17, FontFamily -> fontfamily}; style2 = Directive[Black(*,Thick*)];
	kbdat = Transpose[{ConstantArray[(*kdat*)Range[n], m], bdat}, {3, 1, 2}];
	colors = Map[cfunc, Rescale @ cdat, {2}];
	lines = MapThread[If[OptionValue[Joined], Line, Point][#, VertexColors -> #2] &, {kbdat, colors}];
	(*ps1 = Sequence @@ FilterRules[{ps}, Options[Graphics]]; ps2 = Sequence @@ FilterRules[{ps}, Options[BarLegend]];*)
	ps1 = optionsselect[ps, Graphics];
	ps2 = optionsselect[ps, BarLegend];
	bfig = Graphics[{Thick, lines}, ps1, GridLines -> {ptsnumbers, Automatic}, PlotRangeClipping -> True, (*AspectRatio -> GoldenRatio,*) FrameTicks -> frameticks, 
					 Frame -> True, FrameLabel -> {None, "\!\(\*SubscriptBox[\(E\), \(\[VeryThinSpace]\)]\)"}, FrameTicksStyle -> style2, FrameStyle -> style2, LabelStyle -> style2, BaseStyle -> style];
	legend = BarLegend[{cfunc, {0, 1}}, ps2, Ticks -> Transpose[{{0, 1}, NumberForm[#, {3, 4}] & /@ MinMax[cdat]}], (*"Ticks" -> {0, 1}, "TickLabels" -> {"Min", "Max"},*) TicksStyle -> style2, FrameStyle -> style2, LabelStyle -> style];
	Legended[bfig, legend]
];*)


