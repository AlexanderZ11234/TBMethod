(* ::Package:: *)

BeginPackage["TBMethod`", 
			   {"TBMethod`MDConstruct`",
			    "TBMethod`EigenSpect`",
				"TBMethod`LGFF`",
				"TBMethod`DataVisualization`"}
			]

(* Declare your package's public symbols here. *)

SayHello

Begin["`Private`"]

(* Define your public and private symbols here. *)

$DistributedContexts = {"Global`", "TBMethod`"};

$PlotTheme = {"Scientific", "LargeLabels", "SansLabels"};

SayHello[name_?StringQ] := Print["Hello ", name, "!"]


End[] (* End `Private` *)

EndPackage[]
