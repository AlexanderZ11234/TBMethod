(* ::Package:: *)

(* ::Text:: *)
(*This file collects old functions.*)
(*A function being old here almost means it is created earlier in the D&R phase but works well but not as efficient as the latest version in the runtime files.*)
(*But hidden bugs and errors can not be 100% ruled out, and because they are old, I will not devote more effort to repair them.*)
(*Or, at least the newer versions are their efficiency reparations.*)


(* ::Title:: *)
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




(* ::Title:: *)
(*EigenSpect`**)


(* ::Title:: *)
(*LGFF`**)


(* ::Title:: *)
(*DataVisalization`**)
