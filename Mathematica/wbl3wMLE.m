(* ::Package:: *)

(* ::Title:: *)
(*Fitting the Weibull distribution using weighted MLE*)


(*: ====================================================*)
(*: NAME:   wMLE`                           *)
(*: Context:wMLE`                           *)
(*: ====================================================*)
(*: Author: Denis Cousineau                             *)
(*: Email:  Denis.Cousineau@uottawa.ca                  *)
(*: ====================================================*)
(*: Title: Tools for fitting data using the weighted   *)
(*: Maximum Likelihood estimator.                       *)
(*: Reference: Cousineau, D. (2009). Nearly unbiased estimators *)
(*:     for the three\[Hyphen]parameter Weibull distribution with greater *)
(*:     efficiency than the iterative likelihood method. British Journal *) 
(*:     of Mathematical and Statistical Psychology, 62(1), 167-191.     *)
(*: ====================================================*)
(*: Summary:This package provides a function,           *)
(*  wMLE, which peforms MLE of the 3-parameter *)
(*  Weibull distribution. It uses weights     *)
(*: J1, J2 & J3, but other weights can be specified.    *)
(*: ====================================================*)


BeginPackage["wMLE`"]


(* ::Section:: *)
(*(*: ================================*)*)
(*(*: Usage and error messages                           *)*)
(*(*: ================================*)*)


(* ::Subsection:: *)
(*Usage for the defined functions*)


J1::usage="J1[n] returns the median weight for a given sample size n. This quantity does not depend on the true parameters (it equals the median of a Gamma[n,1/n])."
J2::usage="J2[n] returns the median weight for a given sample size n. This quantity does not depend on the true parameters (obtained through Monte Carlo simulations)."
J3::usage="J3[n,g] returns the median weight for a sample size n and a shape parameter of g (obtained through Monte Carlo simulations)."


SetHeuristics::usage="SetHeuristics[{\!\(\*SubscriptBox[\"\[Gamma]\", \"lo\"]\),\!\(\*SubscriptBox[\"\[Gamma]\", \"hi\"]\)},{\!\(\*SubscriptBox[\"\[Alpha]\", \"lo\"]\),\!\(\*SubscriptBox[\"\[Alpha]\", \"hi\"]\)}] set the starting value intervals for \[Gamma] and \[Alpha] in the search. Defaults are {0.6,2.4},{280,600};
SetHeuristics[max\[Gamma]] set the maximum allowed value for the shape parameter. Default is 5."


wMLE::usage="wMLE[data] performs a ML search using the weights explicitely given. Use SetHeuristics to prepare 
the starting values.\n The weights herein are the J1 to J3 weights (in the median sense).\n\t"


(* ::Subsection:: *)
(*Error messages*)


SetHeuristics::invneg="The lower bound for \[Gamma] cannot be smaller than zero.";
SetHeuristics::invbnds="Invalid bounds. The lower bound `1` must be smaller than the upper `2`";
SetHeuristics::invmax="The maximum allowed value `1` for \[Gamma] must be positive.";

wMLE::nosol="The solution found exceed 0.001; you may be in a region with no consistant solution (\[Gamma] \!\(\*OverscriptBox[\"<\", \"?\"]\) 1).";
wMLE::badbnd="Bad starting value heuristics; use SetHeuristics.";
wMLE::invbnds="Invalid bounds. The lower bound `1` must be smaller than the upper `2`"
wMLE::invneg="The lower bound for \[Gamma] cannot be smaller than zero."
wMLE::invmax="The maximum allowed value `1` for \[Gamma] must be positive."

J1::notcomp="The lookup table for J1 does not include n = `1`. Using asymptotic value 1.";
J2::notcomp="The lookup table for J2 does not include n = `1`. Using asymptotic value 1.";
J3::notcomp="The lookup table for J3 does not include n = `1`. Using asymptotic value \[Gamma]/(\[Gamma]-1).";


(* ::Section:: *)
(*(*: =================================*)*)
(*(*: Lookup tables                                              *)*)
(*(*: =================================*)*)


Begin["`Private`"]


(* ::Subsection:: *)
(*Loading lookup tables*)


(* ::Subsubsection:: *)
(*Lookup tables for J1, J2, J3.*)


AllJ1=Import[NotebookDirectory[]<>"\\..\\weigths\\J1.tsv","TABLE"];
AllJ2=Import[NotebookDirectory[]<>"\\..\\weigths\\J2.tsv","TABLE"];
AllJ3=Import[NotebookDirectory[]<>"\\..\\weigths\\J3.tsv","TABLE"];


(* ::Section:: *)
(*(*: =================================*)*)
(*(*: Code                                                          *)*)
(*(*: =================================*)*)


(* ::Subsection:: *)
(*Various initialization*)


Off[NIntegrate::"inum"];
Off[CompiledFunction::"cfsa"];
Off[NIntegrate::"ncvb"];
Off[NIntegrate::"slwcon"];

nofit={-1,{\[Gamma]->-1,\[Beta]->-1,\[Alpha]->-1}};
m\[CapitalDelta]shape=0.25;p\[CapitalDelta]shape=2.75;
m\[CapitalDelta]shift=280;p\[CapitalDelta]shift=600;
MaxG=5;

AllJ1Known=Union[AllJ1[[All,1]]]
AllJ2Known=Union[AllJ2[[All,1]]]
AllJ3Known=Union[AllJ3[[All,1]]]


(* ::Subsection:: *)
(*Defining the median weights J1, J2 and J3 (two of them thru look-up tables)*)


J1[n_Integer]:=If[MemberQ[AllJ1Known,n],
First[Select[AllJ1,#[[1]]==n&]][[2]],
Message[J1::notcomp,n];1
];
J2[n_Integer]:=If[MemberQ[AllJ2Known,n],
First[Select[AllJ2,#[[1]]==n&]][[2]],
Message[J2::notcomp,n];1
];


J3[n_Integer,g_]:=Module[{sol,myg,ming,maxg, resolution},
	myg=g;
	resolution=AllJ3[[2,2]]-AllJ3[[1,2]];
	ming=Min[AllJ3[[All,2]]];maxg=Max[AllJ3[[All,2]]];
	If[myg>maxg,myg=maxg];
	If[myg<ming,myg=ming];

	If[Not[MemberQ[AllJ3Known,n]],
		Message[J3::notcomp,n];myg/(myg-1),

		sol=Select[AllJ3,#[[1]]==n&&Abs[#[[2]]-myg]<resolution&];
		If[Length[sol]==1,
			First[sol][[3]],
			((myg-sol[[1,2]])/resolution)sol[[2,3]]+((sol[[2,2]]-myg)/resolution)sol[[1,3]]
		]
	]
]


(* ::Subsection:: *)
(*Various procedure*)


SetHeuristics[{m\[CapitalDelta]\[Gamma]_,p\[CapitalDelta]\[Gamma]_},{m\[CapitalDelta]\[Alpha]_,p\[CapitalDelta]\[Alpha]_}]:=Module[{},
	If[m\[CapitalDelta]\[Gamma]>=p\[CapitalDelta]\[Gamma],Message[SetHeuristics::invbnds,m\[CapitalDelta]\[Gamma],p\[CapitalDelta]\[Gamma]]];
	If[m\[CapitalDelta]\[Alpha]>=p\[CapitalDelta]\[Alpha],Message[SetHeuristics::invbnds,m\[CapitalDelta]\[Alpha],p\[CapitalDelta]\[Alpha]]];If[m\[CapitalDelta]\[Gamma]<= 0,Message[SetHeuristics::invneg,m\[CapitalDelta]\[Gamma]]];

	If[m\[CapitalDelta]\[Gamma]<p\[CapitalDelta]\[Gamma]&&m\[CapitalDelta]\[Alpha]<p\[CapitalDelta]\[Alpha]&&m\[CapitalDelta]\[Gamma]>0,
		m\[CapitalDelta]shape=m\[CapitalDelta]\[Gamma];p\[CapitalDelta]shape=p\[CapitalDelta]\[Gamma];
		m\[CapitalDelta]shift=m\[CapitalDelta]\[Alpha];p\[CapitalDelta]shift=p\[CapitalDelta]\[Alpha];
	]
]
SetHeuristics[max\[Gamma]_]:=Module[{},
	If[max\[Gamma]<0,
		Message[SetHeuristics::invmax,max\[Gamma]],
		MaxG=max\[Gamma]
	];
]


(* ::Subsection:: *)
(*The algebraic solutions*)


fctGamma[x_,n_,{\[Gamma]_,\[Alpha]_,w2_}]:=w2/\[Gamma]-(\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(\((
\*SuperscriptBox[\((x[\([i]\)] - \[Alpha])\), \(\[Gamma]\)]*Log[x[\([i]\)] - \[Alpha]])\)/\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]
\*SuperscriptBox[\((x[\([i]\)] - \[Alpha])\), \(\[Gamma]\)]\)\)\))+1/n (\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]\(Log[x[\([i]\)] - \[Alpha]]\)\))


fctBeta[x_,n_,{\[Gamma]_,\[Alpha]_,w1_}]:=Power[1/(n w1) (\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]
\*SuperscriptBox[\((x[\([i]\)] - \[Alpha])\), \(\[Gamma]\)]\)), (\[Gamma])^-1]


fctAlpha[x_,n_,{\[Gamma]_,\[Alpha]_,w3_}]:=(1/n)*(\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]
\*FractionBox[\(1\), \(x[\([i]\)] - \[Alpha]\)]\))*(\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]
\*SuperscriptBox[\((x[\([i]\)] - \[Alpha])\), \(\[Gamma]\)]\))/(\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(n\)]
\*SuperscriptBox[\((x[\([i]\)] - \[Alpha])\), \(\[Gamma] - 1\)]\))-w3


(* ::Subsection:: *)
(*The fitting procedures per se*)


wMLE[sample_List, opt___]:=Module[{ss,sol1,sol2,mil, \[Alpha],\[Gamma],nmrules},

nmrules=FilterRules[{opt},Options[NMinimize]];

ss=Length[sample];
mil=(m\[CapitalDelta]shape+p\[CapitalDelta]shape)/2;

If[Min[sample]<m\[CapitalDelta]shift,
Message[wMLE::badbnd];p\[CapitalDelta]shift=Min[sample];m\[CapitalDelta]shift=0.9Min[sample]
];

Quiet[
(*the objective function*)
obj[pg_?NumericQ,pa_?NumericQ]:=((fctGamma[sample,ss,{pg,pa,J2[ss]}]^2)+(fctAlpha[sample,ss,{pg,pa,J3[ss,pg]}]^2));

sol1=NMinimize[
	{obj[\[Gamma],\[Alpha]], 0 < \[Gamma] < MaxG && -\[Infinity] < \[Alpha] <Min[sample]},
	{{\[Alpha],m\[CapitalDelta]shift,p\[CapitalDelta]shift},{\[Gamma],m\[CapitalDelta]shape,p\[CapitalDelta]shape}},
	Sequence[nmrules],
	Method->{"NelderMead","PostProcess"->False},
	AccuracyGoal->10,
	MaxIterations->200,
	EvaluationMonitor:>(Nothing;
		(*Print[{\[Gamma],\[Alpha]}, (fctGamma[sample,ss,{\[Gamma],\[Alpha],wei2}]^2+fctAlpha[sample,ss,{\[Gamma],\[Alpha],wei3}]^2)];*)
	)
	]
];

sol2=fctBeta[sample,ss,{\[Gamma]/.sol1[[2]],\[Alpha]/.sol1[[2]],J1[ss]}];

If[sol1[[1]]>0.001,Message[wMLE::nosol]];
{sol1[[1]]//Chop,
{Global`\[Gamma]->\[Gamma]/.sol1[[2]],Global`\[Beta]->sol2,Global`\[Alpha]->\[Alpha]/.sol1[[2]]}
}
]


(* ::Section:: *)
(*(*: ==============================*)*)
(*(*: Closing instructions                               *)*)
(*(*: ==============================*)*)


End[]


SetAttributes[{J1,J2,J3},ReadProtected];
SetAttributes[{SetHeuristics},ReadProtected];
SetAttributes[wMLE,ReadProtected];


CellPrint[Cell[BoxData[{
StyleBox[RowBox[{"Package wMLE`"}],FontFamily->"Garamond",FontSize->16],
StyleBox[RowBox[{"\[Copyright] Denis Cousineau (2007), denis.cousineau@uottawa.ca"}],FontFamily->"Garamond"],
StyleBox[RowBox[{"This software is free of charge for academics and students"}],FontFamily->"Garamond",FontWeight->"Plain"],
StyleBox[RowBox[{"Use ? wMLE`* to see the list of functions available."}],FontFamily->"Garamond",FontWeight->"Plain"]
}],
"Output",
CellFrame->{{0,0},{2,2}},Background->GrayLevel[.8]
]]


EndPackage[]


(* ::Section:: *)
(*(*: ================================*)*)
(*(*: End of package                                          *) *)
(*(*: ================================*)*)
