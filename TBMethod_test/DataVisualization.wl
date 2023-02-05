(* ::Package:: *)

BeginPackage["TBMethod`DataVisualization`"]

(* Declare your package's public symbols here. *)

(* Exported symbols added here with SymbolName::usage *)
BandPlot::usage = "Plots band structures with high symmetry points annotated.";

Vector2DPlot::usage = "Visualizes a 2D vector field with local magnitudes and an emphasis on local directions.";

Stream2DPlot::usage = "Visualizes a 2D vector field with local magnitudes and an emphasis on nonlocal tendencies of local directions.";

LocalDOSPlot::usage = "Visualizes local density of states in either real or reciprocal space.";

LocalDOSTidy::usage = "Tidies up local density of states for a better effect of visualization."

Begin["`Private`"]
(* Implementation of the package *)
(*SetOptions[{ParallelSum}, Method -> "ItemsPerEvaluation" -> 100 $KernelCount];*)
$DistributedContexts = {"Global`", "TBMethod`"}; (*Otherwise, dim in HMatrixFromHoppings is not working in parallel subkernels. This is very IMPORTANT!!!*)


(*BandPlot[banddata_, hisymmptname: {_String..}, ptsnumbers: {_Integer..}, s:OptionsPattern[ListLinePlot]] :=*)
BandPlot[banddata_, hisymmptname: {_String...}: {""}, ptsnumbers: {_?NumericQ...}: {}, s:OptionsPattern[ListLinePlot]] :=
Module[{(*plottheme, *)frameticks, gridlines, dticks},
	(*plottheme = {"Scientific", "SansLabels", "LargeLabels"};*)
	dticks = {ptsnumbers, hisymmptname}\[Transpose];
	frameticks = {{Automatic, None}, {dticks, None}};
	ListLinePlot[banddata, s, 
		PlotStyle -> Blue, AspectRatio -> 1.5, (*PlotTheme -> plottheme,*)
		FrameTicks -> frameticks, GridLines -> {ptsnumbers, Automatic}
	]
];

Vector2DPlot[jvekdata_, ops:OptionsPattern[ListVectorDensityPlot]] :=
ListVectorDensityPlot[
	jvekdata, ops, AspectRatio -> Automatic, ColorFunction -> "SunsetColors", MaxRecursion -> 2,
	VectorScaling -> "Log", VectorSizes -> {0, 1}, VectorColorFunction -> None, InterpolationOrder -> 1,
	VectorStyle -> Green(*Black*), VectorPoints -> "Hexagonal", VectorMarkers -> "PinDart"
];
	
Stream2DPlot[jvekdata_, ops:OptionsPattern[ListStreamDensityPlot]] :=
ListStreamDensityPlot[
	jvekdata, ops, AspectRatio -> Automatic, ColorFunction -> "SunsetColors", StreamScale -> {.2, All, Scaled[.15]},
	StreamColorFunction -> None, StreamStyle -> Green, StreamMarkers -> "Dart", StreamPoints -> Fine,
	InterpolationOrder -> 1
];

LocalDOSPlot[data_, ops:OptionsPattern[ListDensityPlot]] :=
ListDensityPlot[
	data, ops, AspectRatio -> 1.5, PlotRange -> All,
	ColorFunction -> (ColorData["SunsetColors"][#^(1/2)]&)
];


LocalDOSTidy[data_, quantile_] :=
Module[{maxquant, clipped, min = Min[data]},
	maxquant = Quantile[data // Flatten, quantile];
	clipped = Clip[data, {min, maxquant}];
	GaussianFilter[clipped, 2]
];


End[] (* End `Private` *)

EndPackage[]
