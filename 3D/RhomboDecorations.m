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

Combine::usage=
"
"

CreateMasterList::usage=
"
"

DecorateGrid::usage =
"DecorateGrid[til_,decorations_]
"

GenerateRotatedTile::usage =
"GenerateRotatedTile[triple]
"

GenerateRotatedTiles::usage =
"GenerateRotatedTiles[triples]
"

GetTriples::usage =
"GetTriples[]
"

PlotTiles::usage =
"PlotTile[tiles]
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
			TransformRhombs,
			Combine,
			CreateMasterList,
			PrintLengths,
			DecorateGrid,
			GenerateRotatedTile,
			GenerateRotatedTiles,
			GetTriples,
			PlotRotatedTiles]

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
		TransformRhombs,
		Combine,
		CreateMasterList,
		PrintLengths,
		DecorateGrid,
		GenerateRotatedTile,
		GenerateRotatedTiles,
		GetTriples,
		PlotRotatedTiles,
		DecorateGrid,
		GenerateRotatedTile,
		GenerateRotatedTiles,
		GetTriples,
		PlotRotatedTiles]

Needs["TetGenLink`"]

$filenameoblate = "oblate-decoration.csv"
$filenameprolate = "prolate-decoration.csv"
$filenamedecoratedrhombs = "rhombs-decoration-j.csv"
$tol=10^-6
$round=10^-9
$triples = Select[Tuples@{Range[6],Range[6],Range[6]},Less@@#&];

$tau = (1+Sqrt[5])/2;
$thta0 = ArcCos[1/(2$tau-1)];(* = ArcCos[1/Sqrt[5]] for QC*)
$thta1 = ArcCos[($tau-1)/2];(* = 2Pi/5 for QC*)
$r = Union[{{0,0,1}},
				Table[CoordinateTransform["Spherical"->"Cartesian",{1,$thta0,x}],
				{x,{-2$thta1,-$thta1,0,$thta1,2$thta1}}]];

  
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

CreateMasterList[decoratedtiles_]:=
Module[{n,dectiles,counter=0,i,masterlist={}},

		n=Length[decoratedtiles];

		dectiles=decoratedtiles;
		
		(*combines all subgraphs into one largest list*)
		(*increments the indices accordingly*)
		
		Do[
			dectiles[[i,All,1]]=Plus[dectiles[[i,All,1]],counter];
			dectiles[[i,All,3]]=Plus[dectiles[[i,All,3]],counter];
			counter=counter+Length[dectiles[[i]]];
			masterlist=Join[masterlist,dectiles[[i]]];,
		{i,n}];

		masterlist=Chop[masterlist];
		masterlist=Round[masterlist,$round];

		(*masterlist=SortBy[SortBy[SortBy[masterlist,#[[2,3]]&],#[[2,2]]&],#[[2,1]]&];*)
		masterlist=SortBy[masterlist,#[[2]]&];

		masterlist=N[masterlist];
		Print[Length[Union[masterlist[[All,2]]]]];

		Print["finished master list"];

	masterlist
]

Combine[decoratedtiles_]:=
Module[{n,counter=0,i,masterlist={},newlist={},keys,vals,dectiles,CheckEqual,
		coords,neighbors,assoc},

		n=Length[decoratedtiles];

		dectiles=decoratedtiles;
		
		masterlist=CreateMasterList[decoratedtiles];		

		n=Length[masterlist];

		CheckEqual[r1_,r2_]:=
			If[Norm[r1-r2]<$tol,
				Return[True];,
				Return[False];
			];

		(*identifies unique nodes, ASSUMING degeneracies are neighboring...*)
		(*this assumption relies on algorithm of Sort[] to work properly*)

		counter=1;
		assoc = Association[masterlist[[1,1]]->1];
		coords={masterlist[[1,2]]};
		neighbors={masterlist[[1,3]]};

		Do[
			If[Mod[i,5000]==0,Print["merging ", i, " of ", n];];
			If[CheckEqual[masterlist[[i-1,2]],masterlist[[i,2]]],
					neighbors[[counter]] = Join[neighbors[[counter]],masterlist[[i,3]]];,
				counter=counter+1;
				AppendTo[coords,masterlist[[i,2]]];
				AppendTo[neighbors,masterlist[[i,3]]];
			];
			AppendTo[assoc,masterlist[[i,1]]->counter];
			,
		{i,2,n,1}];

	neighbors=neighbors/.assoc;
	neighbors=Table[Union[x],{x,neighbors}];
	{coords,neighbors}

]

PrintLengths[decgraph_]:=
Module[{Lengths,curcount,n},
		Lengths=Table[Length[x],{x,decgraph[[2]]}];
		n=Length[Lengths];
		Do[
			curcount=Length[Position[Lengths,i]];
			Print["# vertices with ", i, " neighbor(s) = ", curcount, " (", N[curcount/n]*100,"%)"];,
		{i,Union[Lengths]}
		]
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

DecorateGrid[til_,decorations_]:=
Module[{decoratedtiles,t,add,tripleindex,n,decorationj,T},

			n = Length[til[[5]]];

			add[p_]    := Apply[Plus, p*$r];

			decoratedtiles = Table[tripleindex = Flatten[Position[$triples,til[[5,i]]]][[1]];
									decorationj=Part[decorations,tripleindex];
									T=add[til[[4,i,1]]];
									decorationj[[All,2]] = Table[T+x,{x,decorationj[[All,2]]}];
									decorationj,
								{i,n}];
		decoratedtiles
]

GenerateRotatedTile[triple_]:=
Module[{u,e1,e2,e3,e4,s,c,v,volume,voloblate,volprolate,type,
		t12,t13,t23,T={0,0,0},R1,R2,R3,R2theta,tile,curVertices,vol,n,thta,t,
		pr1={0,0,1},newtile={},i,curStart,curEnd,SkewMatrix,GetTile,AxisAngleMatrix},
		
			Print[triple];

			e1=$r[[triple[[1]]]];
			e2=$r[[triple[[2]]]];
			e3=$r[[triple[[3]]]];

			volume = FullSimplify[Norm[e1.Cross[e2,e3]]];
			voloblate = 1/5 Sqrt[10-2 Sqrt[5]];
			volprolate= 1/5 Sqrt[2 (5+Sqrt[5])];

			If[volume==voloblate,
				type="Oblate";,
				If[volume==volprolate,
					type="Prolate";,type="Unknown";
				];
			];

			SkewMatrix[v_]:=
			Module[{v1=v[[1]],v2=v[[2]],v3=v[[3]]},
					{{0,-v3,v2},
					{v3,0,-v1},
					{-v2,v1,0}}];

			AxisAngleMatrix[n_,thta_]:=
			Module[{x=n[[1]],y=n[[2]],z=n[[3]],c1=Cos[thta],s1=Sin[thta],c2=1-Cos[thta]},
					{{x x c2, x y c2- z s1, x z c2+ y s1},
					{y x c2+ z s1, y y c2, y z c2- x s1},
					{z x c2- y s1, z y c2+ x s, z z c2}}+c1*IdentityMatrix[3]];

			curVertices={{0,0,0},e1,e2,e3,e1+e2,e3+e1,e2+e3,e1+e2+e3};

			GetTile[t_]:=
				{{t[[1]],t[[2]]},{t[[1]],t[[3]]},{t[[2]],t[[5]]},{t[[3]],t[[5]]},{t[[4]],t[[6]]},
				{t[[4]],t[[7]]},{t[[6]],t[[8]]},{t[[7]],t[[8]]},{t[[1]],t[[4]]},{t[[2]],t[[6]]},{t[[3]],t[[7]]},{t[[5]],t[[8]]}};

			tile  = Apply[GetTile,{curVertices}];

			t12=e1.e2;t13=e1.e3;t23=e2.e3;
			If[t12<0&&t13<0,T=-e1;e1=-e1;];
			If[t12<0&&t23<0,T=-e2;e2=-e2;];
			If[t13<0&&t23<0,T=-e3;e3=-e3;];
			If[t12>0&&t23>0,e4=e1;e1=e2;e2=e4;];
			If[t13>0&&t23>0,e4=e1;e1=e3;e3=e4;];

			If[e1.pr1==1,
				R1=IdentityMatrix[3];,
				v=SkewMatrix[Cross[e1,pr1]
			];

			s=Norm[v];
			c=e1.pr1;
			R1=IdentityMatrix[3]+v+v.v(1-c)/s^2];

			If[Cross[e2,e3].e1>0,u=R1.e2,u=R1.e3];

			If[u[[1]]==0,
				If[u[[2]]<0,
					R2theta=Pi/4;
					R2theta=-Pi/4;
				];,
				If[u[[1]]<0,
					R2theta=-ArcTan[u[[2]]/u[[1]]] + Pi;,
					R2theta=-ArcTan[u[[2]]/u[[1]]];
				];
			];

			R2=AxisAngleMatrix[pr1,R2theta];

			R3=R2.R1;

			For[i=1,i<=Length[tile],i++,
				curStart=R3.(T+tile[[i,1]]);
				curEnd=R3.(T+tile[[i,2]]);
				newtile=Append[newtile,{curStart,curEnd}];
			];

	{T,R3,tile,newtile,type}
]

GenerateRotatedTiles[triples_]:=
	Table[Flatten[Map[GenerateRotatedTile,{t}],1],{t,triples}]

GetTriples[]:=$triples

PlotRotatedTiles[tiles_]:=
Module[{oldtiles,newtiles},
		oldtiles = tiles[[All,4]];
		newtiles = tiles[[All,5]];

		Show[Table[Graphics3D[Join[{Black,Map[Line,t]}]],{t,oldtiles}],
			Table[Graphics3D[Join[{Blue,Map[Line,t]}]],{t,newtiles}],
			Graphics3D[Arrow[{{0,0,0},#}]&/@{{1,0,0},{0,1,0},{0,0,2}}, Axes->True],
			Graphics3D[Arrow[{{0,0,0},#}]&/@$r, Axes->True],
		ViewPoint->{5,-1,2}*1000000,
		ImageSize->{700,700}]
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
		TransformRhombs,
		Combine,
		CreateMasterList,
		PrintLengths]

(********************************************************************)
(*                                                                  *)
(*        End of package "Rhombohedra`oblateprolate`"               *)
(*                                                                  *)
(********************************************************************)

EndPackage[]
