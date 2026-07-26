(* ::Package:: *)

(** User Mathematica initialization file **)

(*Block[{kernelnumber, machinenames = {"desktop-p24f3e7"}},
	kernelnumber = If[MemberQ[$MachineName][machinenames], 20, 2 $ProcessorCount];
	LaunchKernels[kernelnumber];
	Echo[StringTemplate["The launched kernel number is: ``."][$KernelCount]];
];*)

LaunchKernels[20];
Echo[StringTemplate["The launched kernel number is: ``."][$KernelCount]];
Needs["TBMethod`"]
ParallelNeeds["TBMethod`"]
Scan[Echo @* Information] @ {"TBMethod`MDConstruct`*", "TBMethod`EigenSpect`*", "TBMethod`LGFF`*", "TBMethod`DataVisualization`*"}
