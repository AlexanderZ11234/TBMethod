(* ::Package:: *)

BeginPackage["TBMethod`EigenSpect`"]

(* Declare your package's public symbols here. *)

(* Exported symbols added here with SymbolName::usage *)
ReciprocalVectors::usage = "Calculates the reciprocal lattice vectors, given the ones in the real space.";

BandData::usage = "Calculates band data of a model over a list of parameter samplings.";

ParallelBandData::usage = "Calculates band data of a model over a list of parameter samplings in parallel manner.";

EigenspectralData::usage = "Calculates eigensystem data of a model over a list of parameter samplings.";

ParallelEigenspectralData::usage = "Calculates eigensystem data of a model over a list of parameter samplings in parallel manner.";

BandDataWithWeight::usage = "Calculates band data with the weight information (e.g., IPR) of a model over a list of parameter samplings.";

ParallelBandDataWithWeight::usage = "Calculates band data with the weight information (e.g., IPR) of a model over a list of parameter samplings in parallel manner.";

PathSample::usage = "Samples evenly points on a path consisting of a sequence of node points, with the the number of points sampled on the shortest segment specified.";

BerryCurvature::usage = "Calculates Berry curvature for an arbitrary matrix model.";

(*FirstBrillouinZone::usage = "Shows the first Brillouin zone with reciprocal lattice vectors.";*)
FirstBrillouinZoneRegion::usage = "Generate the first Brillouin zone as a region.";
FirstBrillouinZonePlot::usage = "Shows the first Brillouin zone with reciprocal lattice vectors.";

LabelPathSamplings::usage = "Labels some specific sampled lattice momentum values, esp. the high-symmetric points.";

(*PlaquetteChern::usage = "Calculates the Chern number contributed by a small region, usually dubbed as a plaquette; summed over the first Brillouin zone cover, the full Chern number is obtained.";
PlaquetteRegionPartition::usage = "Partitions a region, e.g. the first Brillouin zone, with a triangularization cover.";*)

ChernNumberByWilsonLoop::usage = "Calculates the Chern number contributed by a single band or a group of bands, usually the valence bands.";
PlaquetteRegionPartitionComplex::usage = "Partitions a region, e.g. the first Brillouin zone, with a triangularization cover.";

RegionPolarSample::usage = "Samples a 2D region shellwise from its center to its edge.";
ParallelogramCartesianSample::usage = "Samples cartesian-like a paralleogram, which is specificed by its two side vectors.";
WannierChargeCenterByWilsonLoop::usage = "Calculates the Wannier Charge Center contributed by the occupied bands.";

LocalBandGap::usage = "Finds the local band gap among the lattice momentum samplings.";
GlobalBandGap::usage = "Finds the global band gap among the lattice momentum samplings."
ParallelLocalBandGap::usage = "Finds the local band gap among the lattice momentum samplings in parallel.";
ParallelGlobalBandGap::usage = "Finds the global band gap among the lattice momentum samplings in parallel."

BottIndex::usage = "Calculates the Bott index(es)."


Begin["`Private`"]
(* Implementation of the package *)
(*SetOptions[{ParallelSum}, Method -> "ItemsPerEvaluation" -> 100 $KernelCount];*)
$DistributedContexts = {"Global`", "TBMethod`"}; (*Otherwise, dim in HMatrixFromHoppings is not working in parallel subkernels. This is very IMPORTANT!!!*)


ReciprocalVectors[vas_ /; Dimensions[vas] == {2, 2}] :=
Module[{vas3d = BlockDiagonalMatrix[{vas, {{1}}}]},
	Drop[ReciprocalVectors[vas3d], -1, -1]
];
ReciprocalVectors[vas_ /; Dimensions[vas] == {3, 3}] := 2\[Pi] Inverse[vas\[Transpose]];

(*BandData[hbloch_, ks_, map_: Map, s:OptionsPattern[Eigenvalues]] :=
(Sort @ Eigenvalues[hbloch[#], s, Method -> "Direct"] & ~map~ ks)\[Transpose];*)
funcpattern = (_Function | _Symbol | _[__]);

BandData[hbloch:funcpattern, ks_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[hbloch[#], s, Method -> "Direct"] & ~Map~ ks)\[Transpose];
BandData[hbloch:funcpattern, ks_, n_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[hbloch[#], n, s, Method -> "Direct"] & ~Map~ ks)\[Transpose];
ParallelBandData[hbloch:funcpattern, ks_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[hbloch[#], s, Method -> "Direct"] & ~ParallelMap~ ks)\[Transpose];
ParallelBandData[hbloch:funcpattern, ks_, n_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[hbloch[#], n, s, Method -> "Direct"] & ~ParallelMap~ ks)\[Transpose];
(*Version for nonidentity overlapping matrix*)
BandData[{hbloch:funcpattern, sbloch:funcpattern}, ks_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[Through[{hbloch, sbloch}[#]], s, Method -> "Direct"] & ~Map~ ks)\[Transpose];
BandData[{hbloch:funcpattern, sbloch:funcpattern}, ks_, n_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[Through[{hbloch, sbloch}[#]], n, s, Method -> "Direct"] & ~Map~ ks)\[Transpose];
ParallelBandData[{hbloch:funcpattern, sbloch:funcpattern}, ks_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[Through[{hbloch, sbloch}[#]], s, Method -> "Direct"] & ~ParallelMap~ ks)\[Transpose];
ParallelBandData[{hbloch:funcpattern, sbloch:funcpattern}, ks_, n_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[Through[{hbloch, sbloch}[#]], n, s, Method -> "Direct"] & ~ParallelMap~ ks)\[Transpose];


EigenspectralData[hbloch_, ks_, obfunc:funcpattern:Identity, s:OptionsPattern[Eigensystem]] := MapAt[obfunc, {All, 2}][Sort[Eigensystem[hbloch[#], s, Method -> "Direct"]\[Transpose]]] & ~Map~ ks
EigenspectralData[hbloch_, ks_, obfunc:funcpattern:Identity, n_, s:OptionsPattern[Eigensystem]] := MapAt[obfunc, {All, 2}][Sort[Eigensystem[hbloch[#], n, s, Method -> "Direct"]\[Transpose]]] & ~Map~ ks
ParallelEigenspectralData[hbloch_, ks_, obfunc:funcpattern:Identity, s:OptionsPattern[Eigensystem]] := MapAt[obfunc, {All, 2}][Sort[Eigensystem[hbloch[#], s, Method -> "Direct"]\[Transpose]]] & ~ParallelMap~ ks
ParallelEigenspectralData[hbloch_, ks_, obfunc:funcpattern:Identity, n_, s:OptionsPattern[Eigensystem]] := MapAt[obfunc, {All, 2}][Sort[Eigensystem[hbloch[#], n, s, Method -> "Direct"]\[Transpose]]] & ~ParallelMap~ ks

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


ParallelBandDataWithWeight::weightfunc = BandDataWithWeight::weightfunc = "The state function `1` should be an anonymous function or a list of anonymous functions.";

Options[ParallelBandDataWithWeight] = Join[Options[Eigensystem], {"StateFunction" -> (Total[Abs[#]^4] &)}];
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

ParallelBandDataWithWeight[h_, kgrid_, ps:OptionsPattern[]] := ParallelBandDataWithWeight[h, kgrid, Length[h[kgrid[[1]]]], ps]
ParallelBandDataWithWeight[h_, kgrid_, n_Integer, ps:OptionsPattern[]] :=
Module[{ps1val, ps2, eigensystfunc, funcsingle, funcmultiple},
	ps1val = OptionValue["StateFunction"];
	ps2 = TBMethod`DataVisualization`Private`optionsselect[ps, Eigensystem];
	eigensystfunc = k |-> Sort[Eigensystem[h[k], n, ps2, Method -> {"Arnoldi", "MaxIterations" -> \[Infinity]}]\[Transpose]];
	funcsingle = MapAt[ps1val, {All, 2}][eigensystfunc[#]] &;
	funcmultiple = Thread /@ MapAt[Through @* ps1val, {All, 2}][eigensystfunc[#]] &;
	Which[
		MatchQ[ps1val, _Function], ParallelMap[funcsingle, kgrid],
		MatchQ[ps1val, {__Function}], Transpose[ParallelMap[funcmultiple, kgrid], {2, 3, 1}],
		True, Message[ParallelBandDataWithWeight::weightfunc, ps1val]
	]
];

Options[BandDataWithWeight] = Join[Options[Eigensystem], {"StateFunction" -> (Total[Abs[#]^4] &)}];
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
BandDataWithWeight[h_, kgrid_, ps:OptionsPattern[]] := BandDataWithWeight[h, kgrid, Length[h[kgrid[[1]]]], ps]
BandDataWithWeight[h_, kgrid_, n_Integer, ps:OptionsPattern[]] :=
Module[{ps1val, ps2, eigensystfunc, funcsingle, funcmultiple},
	ps1val = OptionValue["StateFunction"];
	ps2 = TBMethod`DataVisualization`Private`optionsselect[ps, Eigensystem];
	eigensystfunc = k |-> Sort[Eigensystem[h[k], n, ps2, Method -> {"Arnoldi", "MaxIterations" -> \[Infinity]}]\[Transpose]];
	funcsingle = MapAt[ps1val, {All, 2}][eigensystfunc[#]] &;
	funcmultiple = Thread /@ MapAt[Through @* ps1val, {All, 2}][eigensystfunc[#]] &;
	Which[
		MatchQ[ps1val, _Function], Map[funcsingle, kgrid],
		MatchQ[ps1val, {__Function}], Transpose[Map[funcmultiple, kgrid], {2, 3, 1}],
		True, Message[bandDataWithWeight::weightfunc, ps1val]
	]
];


PathSample[pts:{{__?NumericQ}..} /; Length[pts] >= 3, nleast_?(# \[Element] PositiveIntegers &)] :=
Module[{ns, normalizeddist, xgridlines, samplings},
	normalizeddist = Normalize[BlockMap[EuclideanDistance @@ # &, pts, 2, 1], Min];
	ns = Floor[normalizeddist nleast];
	xgridlines = Accumulate[Prepend[ns, 1]];
	samplings = Join @@ (Most /@ MapThread[Subdivide, {Most[pts], Rest[pts], ns}]);
	{samplings, xgridlines}
];

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
FirstBrillouinZoneRegion[vbs_List, n_:2] :=
Module[{intcoeffs, reciprocalvectors, len = Length[vbs], longest = Max[Norm /@ vbs], voronoimesh, zero = 1.*^-5},
	intcoeffs = Tuples[Range[-n, n], len];
	voronoimesh=VoronoiMesh[intcoeffs . vbs, ConstantArray[{-1, 1} longest, len]];
	(*MinimalBy[Norm @* RegionCentroid][MeshPrimitives[voronoimesh, len]] // First*)
	SelectFirst[MeshPrimitives[voronoimesh, len], Norm @ RegionCentroid[#] < zero &]
] /; ((Dimensions[vbs] == {2, 2} || Dimensions[vbs] == {3, 3}) && Precision[vbs] == \[Infinity]);

FirstBrillouinZonePlot[vbs_ /; Dimensions[vbs] == {2, 2} || Dimensions[vbs] == {3, 3} , n_:2, opts:OptionsPattern[Show]] :=
Module[{intcoeffs, fbz, reciprocalvectors, len = Length[vbs], \[CapitalGamma], colors},
	intcoeffs = Tuples[Range[-n, n], len]; \[CapitalGamma] = ConstantArray[0, len];
	colors = Take[{Red, Green, Blue}, len]; fbz = FirstBrillouinZoneRegion[vbs, n];
	reciprocalvectors = MapThread[{#, Arrow[{\[CapitalGamma], #2}]} &, {colors, vbs}];
	Show[{HighlightMesh[fbz, {Style[1, Black, Thick], Style[2, Gray, Opacity[0.1]]}],
	If[len == 2, Graphics, Graphics3D][{Thick, reciprocalvectors}]}, opts]
];

(*LabelPathSamplings[pathsamplings_, labels:{__String}] :=
Module[{func, lbllen = Length[labels], numberlen, ptsall, numbers, numbersfinal},
	func = MapAt[lis |-> Callout[lis, #2[[1]], Automatic, Automatic, Appearance -> "CurvedLeader"], #2[[2]]][#] &;
	{ptsall, numbers} = pathsamplings; numberlen = Length[numbers];
	numbersfinal = If[lbllen == numberlen-1, Most @ numbers, MapAt[#-1 &, -1] @ numbers];
	Fold[func, ptsall, {labels, numbersfinal}\[Transpose]]
];*)
LabelPathSamplings[pathsamplings_, labels:{(_String|_OverBar)..}, calloutseqs___, opts:OptionsPattern[Callout]] :=
Module[{func, lbllen = Length[labels], numberlen, ptsall, numbers, numbersfinal},
	func = MapAt[lis |-> Callout[lis, #2[[1]], calloutseqs, opts, Appearance -> "CurvedLeader"], #2[[2]]][#] &;
	{ptsall, numbers} = pathsamplings; numberlen = Length[numbers];
	numbersfinal = If[lbllen == numberlen - 1, Most @ numbers, MapAt[# - 1 &, -1] @ numbers];
	Fold[func, ptsall, {labels, numbersfinal}\[Transpose]]
];


(*Berry's-connection-related topological characteristic method*)
BerryCurvature[veck_, H_, nF_?(# \[Element] PositiveIntegers &), opts:OptionsPattern[Eigensystem]] :=
Module[{vx, vy, eigensyst, \[CapitalOmega]nl, \[Delta]k = 1.*^-7, Hveck = H[veck]},
	{vx, vy} = (Hveck - H[veck - \[Delta]k #])/(*\[HBar]*)\[Delta]k & /@ PauliMatrix[0];
	eigensyst = Sort[Eigensystem[Hveck, opts, Method -> "Direct"]\[Transpose]];
	(*\[CapitalOmega]nl=(#\[LeftDoubleBracket]2\[RightDoubleBracket]\[Conjugate].vx.#2\[LeftDoubleBracket]2\[RightDoubleBracket] #2\[LeftDoubleBracket]2\[RightDoubleBracket]\[Conjugate].vy.#\[LeftDoubleBracket]2\[RightDoubleBracket])/(*\[HBar]^-2*)((#2\[LeftDoubleBracket]1\[RightDoubleBracket]-#\[LeftDoubleBracket]1\[RightDoubleBracket])^2)&;*)(*here in \[LeftDoubleBracket]...\[RightDoubleBracket], "2" for eigenvectors & "1" for eigenvalues*)
	\[CapitalOmega]nl = (#2[[2]]\[Conjugate] . vx . #[[2]] #[[2]]\[Conjugate] . vy . #2[[2]])/(*\[HBar]^-2*)((#2[[1]]-#[[1]])^2)&;(*this \[CapitalOmega]nl has accounted the minus sign automatically*)
	(*-2*)2 Im @ Total[\[CapitalOmega]nl @@@ Tuples[TakeDrop[eigensyst, nF]]]
	(*\[CapitalOmega]nl=(#2\[Conjugate].vx.#4 #4\[Conjugate].vy.#2)/(#3-#)^2&;(*{e1,v1,e2,v2}*)
	-2Total[\[CapitalOmega]nl@@@Join@@@Tuples[TakeDrop[eigensyst,nF]]//Im]*)
];

BerryCurvature[veck_, H_, \[Epsilon]F_?(# \[Element] Reals &), opts:OptionsPattern[Eigensystem]] :=
Module[{vx, vy, eigensyst, \[CapitalOmega]nl, \[Delta]k = 1.*^-7, Hveck = H[veck], group},
	{vx, vy} = (Hveck - H[veck - \[Delta]k #])/(*\[HBar]*)\[Delta]k & /@ PauliMatrix[0];
	eigensyst = Sort[Eigensystem[Hveck, opts, Method -> "Direct"]\[Transpose]];
	(*\[CapitalOmega]nl=(#\[LeftDoubleBracket]2\[RightDoubleBracket]\[Conjugate].vx.#2\[LeftDoubleBracket]2\[RightDoubleBracket] #2\[LeftDoubleBracket]2\[RightDoubleBracket]\[Conjugate].vy.#\[LeftDoubleBracket]2\[RightDoubleBracket])/(*\[HBar]^-2*)((#2\[LeftDoubleBracket]1\[RightDoubleBracket]-#\[LeftDoubleBracket]1\[RightDoubleBracket])^2)&;*)(*here in \[LeftDoubleBracket]...\[RightDoubleBracket], "2" for eigenvectors & "1" for eigenvalues*)
	\[CapitalOmega]nl = (#2[[2]]\[Conjugate] . vx . #[[2]] #[[2]]\[Conjugate] . vy . #2[[2]])/(*\[HBar]^-2*)((#2[[1]]-#[[1]])^2)&;(*this \[CapitalOmega]nl has accounted the minus sign automatically*)
	group = GatherBy[eigensyst, First[#] <= \[Epsilon]F &];
	If[Length[group] == 2, 2 Im @ Total[\[CapitalOmega]nl @@@ Tuples[group]], 0.]
];

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

PlaquetteRegionPartitionComplex[region_, opts:OptionsPattern[TriangulateMesh]] :=
Module[{regiondiscrized, meshcoordinates, plaquettevertexindex, regdim = RegionDimension[region]},
	regiondiscrized = TriangulateMesh[region, opts, 
		MaxCellMeasure -> {regdim -> (RegionMeasure[region].1^regdim)}, Method -> "ConstrainedQuality"];
	meshcoordinates = MeshCoordinates[regiondiscrized];
	(*plaquettevertexindex = List @@@ MeshCells[regiondiscrized, regdim];*)
	plaquettevertexindex = If[regdim == 1, {{DeleteDuplicates @ Flatten[#]}}, #] & [List @@@ MeshCells[regiondiscrized, regdim]];
	Echo[
		MeshRegion[regiondiscrized, PlotTheme -> "Lines", 
		PlotLabel -> StringTemplate["Vertex #: ``, Plaquette #: ``."][Length[meshcoordinates], Length[plaquettevertexindex]]]
	];
	{meshcoordinates, plaquettevertexindex}
];

(*plaquettePhase[occupiedstates_] :=
Module[{stateloop, matD, func},
	stateloop = Partition[occupiedstates, 2, 1, {1, 1}];
	func = {vecs1, vecs2} |-> Outer[#\[Conjugate] . #2 &, vecs1, vecs2, 1];
	matD = Dot @@ func @@@ stateloop;
	Arg @ Det[matD]
];*)

innerProductLoopTensor[occupiedstates_] :=
Module[{stateloop, matD, func},
	stateloop = {#, RotateLeft[#]} & [occupiedstates];
	func = {vecs1, vecs2} |-> Outer[#\[Conjugate] . #2 &, vecs1, vecs2, 1];
	Fold[Dot, MapThread[func][stateloop]] (*Dot @@ has the problem of iteration limit 4096!*)
];

plaquettePhase[occupiedstates_] := Arg @ Det[innerProductLoopTensor[occupiedstates]];
plaquetteBandPhase[occupiedstates_] := Arg @ Eigenvalues[innerProductLoopTensor[occupiedstates]] // Sort;


(*plaquettePhase[occupiedstates_] :=
Module[{stateloop, matD, func},
	stateloop = {#, RotateLeft[#]} &[occupiedstates];
	func = {vecs1, vecs2} |-> Outer[#\[Conjugate] . #2 &, vecs1, vecs2, 1];
	matD = Dot @@ MapThread[func][stateloop];
	Arg @ Det[matD]
];*)

ChernNumberByWilsonLoop[heff_, fbzcomplex:{vks:{__}, inds:{__}}, nF_?(# \[Element] PositiveIntegers &), opts:OptionsPattern[Eigensystem]] :=
Module[{occupiedstatesall, occupiedstatesplaquette},
	occupiedstatesall = Take[Sort[Eigensystem[heff[#], opts, Method -> "Direct"]\[Transpose]], nF][[;;, 2]] & /@ vks;
	(*occupiedstatesall = TakeSmallestBy[Eigensystem[heff[#], opts, Method -> "Direct"]\[Transpose], First, nF][[;;, 2]] & /@ vks;*)
	occupiedstatesplaquette = Extract[occupiedstatesall, inds];
	-(1/(2\[Pi]))Sum[plaquettePhase[occupiedstates], {occupiedstates, occupiedstatesplaquette}]
];

RegionPolarSample[fbzreg_?RegionQ, n_Integer, ptsn_ /; ptsn > 0] :=
Module[{fbzvertices, radialsamplings, azimuthalsample, \[CapitalGamma] = {0., 0.}},
	fbzvertices = SortBy[Arg[Complex @@ #] &] @ PolygonCoordinates[fbzreg];
	radialsamplings = Rest[Subdivide[\[CapitalGamma], #, n]] & /@ fbzvertices;
	azimuthalsample = PathSample[# /. {x_, y___, z_} :> {x, y, z, x}, 1+\[LeftCeiling]ptsn #2[[1]]\[RightCeiling]][[1]] &;
	MapIndexed[azimuthalsample, MapThread[Append, {radialsamplings, fbzvertices}]\[Transpose]]
];

ParallelogramCartesianSample[vbs_List, {m_Integer, n_Integer}] :=
Module[{pts},
	pts=Table[{i, j} . vbs, {i, 0., 1., 1./m}, {j, 0., 1., 1./n}];
	Join[pts, pts[[;;, {1}]], 2]
];

WannierChargeCenterByWilsonLoop[heff_, vksloop_List, nF_?(# \[Element] PositiveIntegers &), opts:OptionsPattern[Eigensystem]] :=
Module[{occupiedstates},
	occupiedstates = Take[Sort[Eigensystem[heff[#], opts, Method -> "Direct"]\[Transpose]], nF][[;;, 2]] & /@ vksloop;
	(*occupiedstates = TakeSmallestBy[Eigensystem[heff[#], opts, Method -> "Direct"]\[Transpose], Re @* First, nF][[;;, 2]] & /@ vksloop;*)
	-1/\[Pi] plaquetteBandPhase[occupiedstates]
];


LocalBandGap[h_, ksamples:{({__Real}|_Real)..}, n_Integer, opts:OptionsPattern[Eigenvalues]] :=
Module[{cveval},
	cveval[k_] := Sort[Eigenvalues[h[k], opts, Method -> "Direct"]][[{n + 1, n}]];
	Min @ Chop @ Table[Subtract @@ cveval[k], {k, ksamples}]
];

GlobalBandGap[h_, ksamples:{({__Real}|_Real)..}, n_Integer, opts:OptionsPattern[Eigenvalues]] :=
Module[{cveval, evals},
	cveval[k_] := Sort[Eigenvalues[h[k], opts, Method -> "Direct"]][[{n + 1, n}]];
	evals = Chop @ Table[cveval[k], {k, ksamples}]\[Transpose];
	Min[evals[[1]]] - Max[evals[[2]]]
];

ParallelLocalBandGap[h_, ksamples:{({__Real}|_Real)..}, n_Integer, opts:OptionsPattern[Eigenvalues]] :=
Module[{cveval},
	cveval[k_] := Sort[Eigenvalues[h[k], opts, Method -> "Direct"]][[{n + 1, n}]];
	Min @ Chop @ ParallelTable[Subtract @@ cveval[k], {k, ksamples}]
];

ParallelGlobalBandGap[h_, ksamples:{({__Real}|_Real)..}, n_Integer, opts:OptionsPattern[Eigenvalues]] :=
Module[{cveval, evals},
	cveval[k_] := Sort[Eigenvalues[h[k], opts, Method -> "Direct"]][[{n + 1, n}]];
	evals = Chop @ ParallelTable[cveval[k], {k, ksamples}]\[Transpose];
	Min[evals[[1]]] - Max[evals[[2]]]
];


BottIndex[unitvasbasis_, pts_, chiralop_, fwithinnerdof:((_Function| _Symbol) | {(_Function|_Symbol)..})][h00_] :=
Module[{uC, H, len = Length[pts], h, q, m, ms, coords = LinearSolve[unitvasbasis\[Transpose], #] & /@ pts, bott},
	uC = Sort[Eigensystem[chiralop, Method -> "Direct"]\[Transpose]][[;;, 2]]; (*eigensystem of the chiral operator*)
	H = SparseArray[uC . h00 . (uC\[ConjugateTranspose])] // Chop;
	h = Partition[H, {1, 1} (2len)][[1, 2]];
	q = # . (#3\[ConjugateTranspose]) & @@ SingularValueDecomposition[h(*,MatrixRank[h],TargetStructure->"Structured"*)](*"compact" SVD of h*);
	bott[m_, q_] := 1/(2\[Pi] I) Tr[MatrixLog[m . q . m\[ConjugateTranspose] . q\[ConjugateTranspose], Method -> "Schur"]];
	(*bott[m_,q_]:=1/(2\[Pi] \[ImaginaryI])Total[Log[Eigenvalues[m.q.m\[ConjugateTranspose].q\[ConjugateTranspose],Method->"Direct"]]];*)
	(*m=KroneckerProduct[DiagonalMatrix[Exp[2\[Pi] \[ImaginaryI] f/@coords]],PauliMatrix[0+3]];*)(*This is wrong!*)
	If[
		Head[fwithinnerdof] === List,
		ms = Table[BlockDiagonalMatrix[ParallelTable[MatrixExp[2\[Pi] I f[coord]], {coord, coords}]], {f, fwithinnerdof}];
		(*ms=BlockDiagonalMatrix/@ParallelTable[MatrixExp[2\[Pi] \[ImaginaryI] f[coord]],{f,fwithinnerdof},{coord,coords}];*)
		bott[#, q] & /@ ms,
		m = BlockDiagonalMatrix[Table[MatrixExp[2\[Pi] I fwithinnerdof[coord]], {coord, coords}]];
		bott[m, q]
	]
];


End[] (* End `Private` *)

EndPackage[]
