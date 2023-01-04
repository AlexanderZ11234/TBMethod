(* ::Package:: *)

BeginPackage["TBMethod`EigenSpect`"]

(* Declare your package's public symbols here. *)

(* Exported symbols added here with SymbolName::usage *)
ReciprocalVectors::usage = "Calculates the reciprocal lattice vectors, given the ones in the real space.";

BandData::usage = "Calculates band data of a model over a list of parameter samplings.";

ParallelBandData::usage = "Calculates band data of a model over a list of parameter samplings in parallel manner.";

EigenspectralData::usage = "Calculates eigensystem data of a model over a list of parameter samplings.";

ParallelEigenspectralData::usage = "Calculates eigensystem data of a model over a list of parameter samplings in parallel manner.";

PathSample::usage = "Samples evenly points on a path consisting of a sequence of node points, with the the number of points sampled on the shortest segment specified.";

BerryCurvature::usage = "Calculates Berry curvature for an arbitrary matrix model.";

(*FirstBrillouinZone::usage = "Shows the first Brillouin zone with reciprocal lattice vectors.";*)
FirstBrillouinZoneRegion::usage = "Generate the first Brillouin zone as a region.";
FirstBrillouinZonePlot::usage = "Shows the first Brillouin zone with reciprocal lattice vectors.";

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
BandData[hbloch_, ks_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[hbloch[#], s, Method -> "Direct"] & ~Map~ ks)\[Transpose];
BandData[hbloch_, ks_, n_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[hbloch[#], n, s, Method -> "Direct"] & ~Map~ ks)\[Transpose];
ParallelBandData[hbloch_, ks_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[hbloch[#], s, Method -> "Direct"] & ~ParallelMap~ ks)\[Transpose];
ParallelBandData[hbloch_, ks_, n_, s:OptionsPattern[Eigenvalues]] := (Sort @ Eigenvalues[hbloch[#], n, s, Method -> "Direct"] & ~ParallelMap~ ks)\[Transpose];

EigenspectralData[hbloch_, ks_, obfunc:(_Function | _Symbol):Identity, s:OptionsPattern[Eigensystem]] := MapAt[obfunc, {All, 2}][Sort[Eigensystem[hbloch[#], s, Method -> "Direct"]\[Transpose]]] & ~Map~ ks
EigenspectralData[hbloch_, ks_, obfunc:(_Function | _Symbol):Identity, n_, s:OptionsPattern[Eigensystem]] := MapAt[obfunc, {All, 2}][Sort[Eigensystem[hbloch[#], n, s, Method -> "Direct"]\[Transpose]]] & ~Map~ ks
ParallelEigenspectralData[hbloch_, ks_, obfunc:(_Function | _Symbol):Identity, s:OptionsPattern[Eigensystem]] := MapAt[obfunc, {All, 2}][Sort[Eigensystem[hbloch[#], s, Method -> "Direct"]\[Transpose]]] & ~ParallelMap~ ks
ParallelEigenspectralData[hbloch_, ks_, obfunc:(_Function | _Symbol):Identity, n_, s:OptionsPattern[Eigensystem]] := MapAt[obfunc, {All, 2}][Sort[Eigensystem[hbloch[#], n, s, Method -> "Direct"]\[Transpose]]] & ~ParallelMap~ ks


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

FirstBrillouinZoneRegion[vbs_ /; Dimensions[vbs] == {2, 2} || Dimensions[vbs] == {3, 3} , n_:2] :=
Module[{intcoeffs, reciprocalvectors, len = Length[vbs]},
	intcoeffs = Tuples[Range[-n, n], len];
	MinimalBy[Norm @* RegionCentroid][MeshPrimitives[VoronoiMesh[intcoeffs . vbs], len]] // First
];

FirstBrillouinZonePlot[vbs_ /; Dimensions[vbs] == {2, 2} || Dimensions[vbs] == {3, 3} , n_:2, opts:OptionsPattern[Show]] :=
Module[{intcoeffs, fbz, reciprocalvectors, len = Length[vbs], \[CapitalGamma], colors},
	intcoeffs = Tuples[Range[-n, n], len]; \[CapitalGamma] = ConstantArray[0, len];
	colors = Take[{Red, Green, Blue}, len]; fbz = FirstBrillouinZoneRegion[vbs, n];
	reciprocalvectors = MapThread[{#, Arrow[{\[CapitalGamma], #2}]} &, {colors, vbs}];
	Show[{HighlightMesh[fbz, {Style[1, Black, Thick], Style[2, Gray, Opacity[0.1]]}],
	If[len == 2, Graphics, Graphics3D][{Thick, reciprocalvectors}]}, opts]
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

(*WannerChargeCenter[] :=.*)




End[] (* End `Private` *)

EndPackage[]
