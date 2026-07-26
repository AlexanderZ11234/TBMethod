(* ::Package:: *)

(** User Mathematica initialization file **)

LaunchKernels[20];
Echo[StringTemplate["The launched kernel number is: ``."][$KernelCount]];
Needs["TBMethod`"]
ParallelNeeds["TBMethod`"]
Scan[Echo @* Information] @ {"TBMethod`MDConstruct`*", "TBMethod`EigenSpect`*", "TBMethod`LGFF`*", "TBMethod`DataVisualization`*"}
