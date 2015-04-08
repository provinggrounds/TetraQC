(* ::Package:: *)

(********************************************************************)
(*                                                                  *)
(* :Title:       Oblate                                             *)
(*                                                                  *)
(* :Authors:     Chaney Lin                                         *)
(*                                                                  *)
(* :Date:        March 31, 2015                                     *)
(*                                                                  *)
(* :Summary:     Creates the oblate decoration                      *)
(*                                                                  *)
(********************************************************************)

BeginPackage["Rhombohedra`oblateprolate`"]

(********************************************************************)
(*                                                                  *)
(*              Usage messages for exported functions               *)
(*                                                                  *)
(********************************************************************)

MakeLines::usage =
"MakeLines[list]
"

CombineLists::usage=
"CombineLists[listA,listB]
"

ReadOblate::usage=
"ReadOblate[]
"

ReadProlate::usage=
"ReadProlate[]
"

ReadDecorated::usage=
"ReadDecorated[]
"

CombineDecorations::usage=
"CombineDecorations[decoratedtiles]
"

ShowRhomb::usage=
"ShowRhomb[rhomb]
"

TransformRhomb::usage=
"TransformRhomb[T,R,rhomb] returns, for every point r in Rhomb, T+R.r
"

TransformRhombs::usage=
"TransformRhombs[listT,listR,rhomb]
"

(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[MakeLines,
			CombineLists,
			ReadOblate,
			ReadProlate,
			ReadDecorated,
			CombineDecorations,
			ShowRhomb,
			TransformRhomb,
			TransformRhombs]

(********************************************************************)
(*                                                                  *)
(*                     Start of Private context                     *)
(*                                                                  *)
(********************************************************************)

Begin["`Private`"]

(*********************************)
(* Clearing previous definitions *)
(*********************************)

Clear[MakeLines,
		CombineLists,
		ReadOblate,
		ReadProlate,
		ReadDecorated,
		CombineDecorations,
		ShowRhomb,
		TransformRhomb,
		TransformRhombs]

Needs["TetGenLink`"]

$filenameoblate = "oblate-decoration.csv"
$filenameprolate = "prolate-decoration.csv"
$filenamedecoratedrhombs = "rhombs-decoration-j.csv";
  
MakeLines[list_]:=
Module[{Edges={}},
		Do[
			Do[
				AppendTo[Edges,{list[[i,2]],list[[j,2]]}],
			{j,list[[i,3]]}];,
		{i,Length[list]}];
	Edges
]

CombineLists[listA_,listB_]:=
Module[{list1,list2,vertexmapping={},counter,pos,
		curindex,newindex,neighbors,nbrindex},
		list1=listA[[All,{1,2,3}]];
		list2=listB[[All,{1,2,3}]];

		(*produce mapping of vertices in list2 to existing vertices in list1, if overlap*)
		
		counter=Length[list1]+1;(*counter for new indices*)
		Do[
			pos=Position[list1[[All,2]],Part[list2[[i]],2]];
			If[Length[pos]>0,
				AppendTo[vertexmapping,{i,list1[[First[First[pos]],1]],list2[[i,3]]}];,(*gives existing index*)
				AppendTo[vertexmapping,{i,counter,list2[[i,3]]}];(*gives new index*)
				counter=counter+1;
			];,
		{i,Length[list2]}];

		(*apply map*)
		Do[
			curindex = i[[1]];
			newindex=i[[2]];
			neighbors=i[[3]];

			(*replaces newindex of current cell*)
			list2 = ReplacePart[list2,{curindex,1}->newindex];

			(*replaces newindex for all neighbors of current cell*)
			Do[
				nbrindex=First[First[Position[list2[[j,3]],curindex]]];
				list2[[j,3]]=Sort[ReplacePart[list2[[j,3]],{nbrindex->newindex}]];,
			{j,neighbors}],
		{i,vertexmapping}];

		(*merge lists*)
		Do[
			curindex=list2[[i,1]];
			pos=Position[list1[[All,1]],curindex];
			If[Length[pos]>0,
				pos=First[First[pos]];
					list1[[pos,3]]=Union[list1[[pos,3]],list2[[i,3]]];,
				AppendTo[list1,list2[[i]]];
			];,
		{i,Length[list2]}];

	list1
]

CombineDecorations[decoratedtiles_]:=
Module[{n,combinedtiles={}},

		n=Length[decoratedtiles];
	
		combinedtiles=decoratedtiles[[1]];

		Do[
			Print["adding tile ", i, " of ", n];
			combinedtiles=CombineLists[combinedtiles,decoratedtiles[[i]]];,
		{i,n}];

	combinedtiles
]


ReadOblate[]:= ToExpression[Import[$filenameoblate,"CSV"]]

ReadProlate[]:= ToExpression[Import[$filenameprolate,"CSV"]]

ReadDecorated[]:=ToExpression[Import[$filenamedecoratedrhombs,"CSV"]]

ShowRhomb[rhomb_]:=
Show[Graphics3D[{Arrow[{{0,0,0},{0,0,2.5}}],Arrow[{{0,0,0},{1,0,0}}],Arrow[{{0,0,0},{0,1,0}}],
	Line[MakeLines[rhomb]],
	{PointSize[0.01],Blue,Point[rhomb[[Flatten[Position[rhomb[[All,4]],4]],2]]]},
	{PointSize[Large],Red,Point[rhomb[[Flatten[Position[rhomb[[All,4]],3]],2]]]},
	{PointSize[Large],Green,Point[rhomb[[Flatten[Position[rhomb[[All,4]],2]],2]]]}
	}],
	ViewPoint->{5,-1,2}*10000000,
	ImageSize->{700,700}
]

TransformRhomb[T_,R_,rhomb_]:=
Module[{IR,newrhomb},
		IR=Inverse[R];
		newrhomb=rhomb;
		newrhomb[[All,2]]=Table[-T+IR.r,{r,rhomb[[All,2]]}];
	newrhomb
]

TransformRhombs[RotatedTiles_,Oblate_,Prolate_]:=
Module[{n,T,R,type,i,rhombs={}},
		n = Length[RotatedTiles];
		Do[
			Print["transforming ", i, " of ", n];
			T=RotatedTiles[[i,1]];
			R=RotatedTiles[[i,2]];
			type=RotatedTiles[[i,5]];
			If[type=="Prolate",
				rhombs = Append[rhombs,TransformRhomb[T,R,Prolate]];
			]
			If[type=="Oblate",
				rhombs = Append[rhombs,TransformRhomb[T,R,Oblate]];
			],
		{i,n}];
	rhombs
]


(********************************************************************)
(*                                                                  *)
(*                      End of Private context                      *)
(*                                                                  *)
(********************************************************************)

End[]

(********************************************************************)
(*                                                                  *)
(*                    Protect exported functions                    *)
(*                                                                  *)
(********************************************************************)

Protect[MakeLines,
		CombineLists,
		ReadOblate,
		ReadProlate,
		ReadDecorated,
		CombineDecorations,
		ShowRhomb,
		TransformRhomb,
		TransformRhombs]

(********************************************************************)
(*                                                                  *)
(*        End of package "AperiodicTilings`GridMethod`"             *)
(*                                                                  *)
(********************************************************************)

EndPackage[]
