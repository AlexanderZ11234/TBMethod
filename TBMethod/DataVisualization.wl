(* ::Package:: *)

BeginPackage["TBMethod`DataVisualization`"]

(* Declare your package's public symbols here. *)

(* Exported symbols added here with SymbolName::usage *)
BandPlot::usage = "Plots band structures with high symmetry points annotated.";

Vector2DPlot::usage = "Visualizes a 2D vector field with local magnitudes and an emphasis on local directions.";

Stream2DPlot::usage = "Visualizes a 2D vector field with local magnitudes and an emphasis on nonlocal tendencies of local directions.";

LocalDOSPlot::usage = "Visualizes local density of states in either real or reciprocal space.";

LocalDOSTidy::usage = "Tidies up local density of states for a better effect of visualization.";

BandPlotWithWeight::usage = "Plots band structures with high symmetry points annotated and with extra weight information obtainable from the corresponding eigenvectors, for example, IPR.";

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

Options[BandPlotWithWeight] = Join[Options[Graphics], Options[BarLegend], {Joined -> True, ColorFunction -> (Hue[2(1 - #)/3] &)}];
(*bandPlotWithWeight[banddatawithstate_,cfunc_,cname_String,joined_:(True|False),ps:OptionsPattern[Graphics]]:=*)
BandPlotWithWeight[banddatawithweight_, hisymmptname: {_String...}: {""}, ptsnumbers: {_?NumericQ...}: {}, ps:OptionsPattern[]] :=
Module[{kbdat, colors, m, n, lines, bfig, legend, fontfamily = (*"Helvetica"*)(*"Times New Roman"*)"Arial", style,
		style2, bdat, cdat, cfunc = OptionValue[ColorFunction], dticks, frameticks, ps1, ps2},
	{bdat, cdat} = Transpose[banddatawithweight, {3, 2, 1}]; {m, n} = Dimensions[bdat];
	dticks = {ptsnumbers, hisymmptname}\[Transpose]; frameticks = {{Automatic, None}, {dticks, None}};
	style = {FontSize -> 17, FontFamily -> fontfamily}; style2 = Directive[(*Black*)Gray(*,Thick*)];
	kbdat = Transpose[{ConstantArray[(*kdat*)Range[n], m], bdat}, {3, 1, 2}];
	colors = Map[cfunc, Rescale @ cdat, {2}];
	lines = MapThread[If[OptionValue[Joined], Line, Point][#, VertexColors -> #2] &, {kbdat, colors}];
	ps1 = Sequence @@ FilterRules[{ps}, Options[Graphics]]; ps2 = Sequence @@ FilterRules[{ps}, Options[BarLegend]];
	bfig = Graphics[{Thick, lines}, ps1, GridLines -> {ptsnumbers, Automatic}, PlotRangeClipping -> True,AspectRatio -> GoldenRatio, FrameTicks -> frameticks, 
					 Frame -> True, FrameLabel -> {None, "\!\(\*SubscriptBox[\(E\), \(\[VeryThinSpace]\)]\)"}, FrameTicksStyle -> style2, FrameStyle -> style2, LabelStyle -> style2, BaseStyle -> style];
	legend = BarLegend[{cfunc, {0, 1}},(*4,*)ps2, (*Ticks->Transpose[{{0,1},MinMax[cdat]}],*) "Ticks" -> {0, 1}, "TickLabels" -> {"Min", "Max"}, TicksStyle -> style2, FrameStyle -> style2, LabelStyle -> style];
	Legended[bfig, legend]
];


End[] (* End `Private` *)

EndPackage[]
