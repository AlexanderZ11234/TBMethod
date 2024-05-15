(* ::Package:: *)

BeginPackage["TBMethod`DataVisualization`"]

(* Declare your package's public symbols here. *)

(* Exported symbols added here with SymbolName::usage *)
BandPlot::usage = "Plots band structures with high symmetry points annotated.";

Vector2DPlot::usage = "Visualizes a 2D vector field with local magnitudes and an emphasis on local directions.";

Stream2DPlot::usage = "Visualizes a 2D vector field with local magnitudes and an emphasis on nonlocal tendencies of local directions.";

RealSpaceLocalDOSPlot::usage = "Visualizes local density of states in the real space in the buble style.";

LocalDOSPlot::usage = "Visualizes local density of states in either the real or the reciprocal space continuously by interpolation.";

LocalDOSTidy::usage = "Tidies up local density of states for a better effect of visualization, by smoothing the data via Gaussian filter.";

BandPlotWithWeight::usage = "Plots band structures with high symmetry points annotated and with extra weight information obtainable from the corresponding eigenvectors, for example, IPR.";

AtomPerLayerPlot::usage = "Plots the atom (site) number per layer after the adaptive partitioning of the central scattering region.";

ExportAdaptivePartitionAnimation::usage = "Exports the animation, e.g. a *.gif file, to demonostrate the process of adaptive partition of the CSR.";

CrystalStructurePlot::usage = "Visualizes the crystal structure with Wigner-Seitz primitive cells and lattice vectors.";

Begin["`Private`"]
(* Implementation of the package *)
(*SetOptions[{ParallelSum}, Method -> "ItemsPerEvaluation" -> 100 $KernelCount];*)
$DistributedContexts = {"Global`", "TBMethod`"}; (*Otherwise, dim in HMatrixFromHoppings is not working in parallel subkernels. This is very IMPORTANT!!!*)


(*BandPlot[banddata_, hisymmptname: {_String..}, ptsnumbers: {_Integer..}, s:OptionsPattern[ListLinePlot]] :=*)
BandPlot[banddata_,
		 hisymmptname: {(_String|OverBar[_String])..}: {""},
		 ptsnumbers: {_?NumericQ..}: {1},
		 yticks :{{_, _}..} : Automatic,
		 s:OptionsPattern[ListLinePlot]] :=
Module[{(*plottheme, *)frameticks, gridlines, dticks},
	(*plottheme = {"Scientific", "SansLabels", "LargeLabels"};*)
	dticks = {ptsnumbers, hisymmptname}\[Transpose];
	frameticks = {{(*Automatic*)yticks, None}, {dticks, None}};
	ListLinePlot[banddata, s, FrameStyle -> Directive[Black],
		PlotStyle -> Blue, (*AspectRatio -> GoldenRatio,*) (*PlotTheme -> plottheme,*)
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

RealSpaceLocalDOSPlot[evalandevec_List, ptsdisk:{{_, _, _}..}|{{_, _}..}, region_?RegionQ, innerdof_Integer:1, ops:OptionsPattern[Graphics]] :=
Module[{op, largestcomps, \[Eta] = 1.*^-4, len = Length[ptsdisk], ratio = 2},
	op = If[MatchQ[ptsdisk, {{_, _, _}..}], KeyValueMap[Append] @* (data |-> GroupBy[data, (#[[;;2]] &) -> Last, Total]), Identity];
	largestcomps = Select[Last[#] > \[Eta] &] @ Join[ptsdisk, {BlockMap[Total, Abs[evalandevec[[2]]]^2, innerdof]}\[Transpose], 2];
	Graphics[
	{{Opacity[.1], Green, region},
	 {Opacity[.4], Red, Disk[{#, #2}, Sqrt[len]ratio #3] & @@@ op[largestcomps]}},
	ops, PlotLabel -> StringTemplate["\!\(\*SubscriptBox[\(E\), \(\[VeryThinSpace]\)]\) = ``"][evalandevec[[1]]]]
] /; (Length[Partition[evalandevec[[2]], innerdof]] == Length[ptsdisk]);

LocalDOSPlot[data_, ops:OptionsPattern[ListDensityPlot]] :=
ListDensityPlot[
	data, ops, PlotRange -> All, FrameStyle -> Directive[Black],
	ColorFunction -> (ColorData["SunsetColors"][#^(1/2)] &)
];

LocalDOSPlot[data_, hisymmptname: {(_String|OverBar[_String])...}: {""}, ptsnumbers: {_?NumericQ...}: {}, ops:OptionsPattern[ListDensityPlot]] :=
Module[{frameticks, gridlines, dticks},
	dticks = {ptsnumbers, hisymmptname}\[Transpose];
	frameticks = {{Automatic, None}, {dticks, None}};
	ListDensityPlot[
		data, ops, PlotRange -> All, FrameStyle -> Directive[Black],
		FrameTicks -> frameticks, GridLines -> {ptsnumbers, Automatic},
		ColorFunction -> (ColorData["SunsetColors"][#^(1/2)] &)
	]
];


LocalDOSTidy[data_, quantile_] :=
Module[{maxquant, clipped, min = Min[data]},
	maxquant = Quantile[data // Flatten, quantile];
	clipped = Clip[data, {min, maxquant}];
	GaussianFilter[clipped, 2]
];

optionsselect[options:Sequence[___Rule]][func_Symbol] := optionsselect[options, func];
optionsselect[options:Sequence[___Rule], func_Symbol] := Sequence @@ FilterRules[{options}, Options[func]];

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

Options[BandPlotWithWeight] = Join[Options[Graphics], Options[BarLegend],
	{Joined -> True, ColorFunction -> (Hue[2(1 - #)/3] &), PlotStyle -> Sequence[Thick, PointSize[.02]], "LegendTickDigits" -> {3, 4}, "LegendPosition" -> Right}];
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
	bfig = Graphics[{OptionValue[PlotStyle], lines}, ps1, GridLines -> {ptsnumbers, Automatic}, PlotRangeClipping -> True, (*AspectRatio -> GoldenRatio,*) FrameTicks -> frameticks, 
					 Frame -> True, FrameLabel -> {None, "\!\(\*SubscriptBox[\(E\), \(\[VeryThinSpace]\)]\)"}, FrameTicksStyle -> style2, FrameStyle -> style2, LabelStyle -> style2, BaseStyle -> style];
	legend = BarLegend[{cfunc, {0, 1}}, ps2, Ticks -> Transpose[{{0, 1}, NumberForm[#, OptionValue["LegendTickDigits"]] & /@ MinMax[cdat]}],
					   (*"Ticks" -> {0, 1}, "TickLabels" -> {"Min", "Max"},*) TicksStyle -> style2, FrameStyle -> style2, LabelStyle -> style];
	Legended[bfig, Placed[legend, OptionValue["LegendPosition"]]]
];

AtomPerLayerPlot[ptscsr_, ptscsraped_Association, ops:OptionsPattern[ListPlot]] :=
Module[{atomnumberperlayer = Length /@ ptscsraped, atomnumbermaintained, totalnum = Length[ptscsr], meannum, stddev, meanweight},
	atomnumbermaintained = Length[ptscsr] == Total[atomnumberperlayer];
	meannum = Mean[atomnumberperlayer] // N; stddev = StandardDeviation[atomnumberperlayer] // N;
	meanweight = ToString[Mean[atomnumberperlayer^3 // N] // ScientificForm, TraditionalForm];
	ListPlot[atomnumberperlayer, ops,
		PlotMarkers -> {"\[FilledCircle]", 10}, Filling -> 0, GridLines -> Automatic, PlotRangePadding -> {Scaled[.05], Scaled[.06]},
		PlotLabel -> StringTemplate["Total Atom #: `` (``)\n\!\(\*OverscriptBox[\(N\), \(_\)]\) = ``, \[Sigma] = ``, \!\(\*OverscriptBox[\(w\), \(_\)]\) = ``"][totalnum, atomnumbermaintained, meannum, stddev, meanweight],
		PlotRange -> {0, Automatic}, FrameLabel -> {"Layer #", "Atom #"}
	]
];

Options[ExportAdaptivePartitionAnimation] = Join[{"DisplayDurations" -> 0.5, "PointSize" -> 0.008}, Options[Graphics]];
ExportAdaptivePartitionAnimation[ptscsr_, ptscsraped_, destpath_String, opts:OptionsPattern[]] :=
Module[{bdrange = (List @@ BoundingRegion[ptscsr])\[Transpose], len = Length[ptscsraped], plot, reversedfigs, figframes, optsgraphics},
	optsgraphics = Sequence @@ FilterRules[{opts}, Options[Graphics]];
	(*plot[pts_]:=ListPlot[pts,opts,PlotMarkers->{{"\[FilledCircle]",7}},PlotStyle->RandomColor[],PlotRange->bdrange,PlotRangePadding->Scaled[.05],ImageSize->Large,AspectRatio->Automatic];*)
	plot[pts_] := Graphics[
		{PointSize[OptionValue["PointSize"]], RandomColor[], Point[pts]}, optsgraphics,
		AspectRatio -> Automatic, PlotRange -> bdrange, PlotRangePadding -> Scaled[.05], ImageSize -> Large];
	reversedfigs = plot /@ Reverse[Values[ptscsraped]];
	figframes = Table[Show[reversedfigs[[ ;; i]]], {i, Range[len]}];
	Export[destpath, figframes, "DisplayDurations" -> OptionValue["DisplayDurations"]]
];

Options[CrystalStructurePlot] = Join @@ (Options /@ {ListPlot, ListPointPlot3D, Graphics, Graphics3D, Show});
CrystalStructurePlot[
	vasbasis_ /; Dimensions[vasbasis] == {2, 2} || Dimensions[vasbasis] == {3, 3},
	crystalstructure:KeyValuePattern[_List -> {__List}] | KeyValuePattern[_List -> {Rule[_, _List]..}],
	n:_Integer:2,
	opts:OptionsPattern[]]:=
Module[{ptsdelablized, pts, cells, latticevectors, len = Length[vasbasis], opts1, opts2, opts3, myListPlot, myGraphics},
	myListPlot = If[len == 2, ListPlot, ListPointPlot3D];
	myGraphics = If[len == 2, Graphics, Graphics3D];
	ptsdelablized = If[MatchQ[crystalstructure, KeyValuePattern[_List -> {Rule[_, _List]..}]], Values /@ crystalstructure, crystalstructure];
	opts1 = (*TBMethod`DataVisualization`Private`*)optionsselect[opts, myListPlot];
	opts2 = (*TBMethod`DataVisualization`Private`*)optionsselect[opts, myGraphics];
	opts3 = (*TBMethod`DataVisualization`Private`*)optionsselect[opts, Show];
	pts = myListPlot[FullSimplify @* KeyMap[LinearSolve[vasbasis\[Transpose]]][ptsdelablized], opts1];
	cells = VoronoiMesh[Tuples[Range@@(n{-1, 1}), len] . vasbasis, PlotTheme -> "Lines", MeshCellHighlight -> {{2, All} -> Opacity[.1], {1, All} -> Blue}];
	latticevectors = myGraphics[{Thick, MapThread[{#, Arrow[{ConstantArray[0, len], #2}]} &, {Take[{Red, Green, Blue}, len], vasbasis}]}, opts2];
	Show[{pts, cells, latticevectors}, opts3, PlotRangeClipping -> True]
];


End[] (* End `Private` *)

EndPackage[]
