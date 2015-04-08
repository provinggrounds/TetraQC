(* ::Package:: *)

(********************************************************************)
(*                                                                  *)
(* :Title:       GridMethod3D                                       *)
(*                                                                  *)
(* :Authors:     Chaney Lin                                         *)
(*                                                                  *)
(* :Date:        April 2015                                     *)
(*                                                                  *)
(* :Summary:     unofficial update generalizes to 3D                *)
(*                                                                  *)
(********************************************************************)

BeginPackage["AperiodicTilings`GridMethod`"]

(********************************************************************)
(*                                                                  *)
(*              Usage messages for exported functions               *)
(*                                                                  *)
(********************************************************************)

PlotGridVectors::usage =
"PlotGridVectors[] plots the grid vectors that generate the icosahedral
quasicrystal"

DualizeGrid::usage =
"DualizeGrid[kmin,kmax,gamma] calculates the icosahedral quasicrystal
using the dual grid method. The initial grid consists of parallel planes
in each direction at positions kmin to kmax. Here, kmin and kmax may either
be given as a list of integers of length n for odd values of n, and n/2 for
even values of n, or as a single integer in which case it will be used for
all grid directions. The optional argument gamma determines the shift, it
can be given in three different ways: either explicitly as a real vector of 
length n for odd values of n, and n/2 for even values of n, or it is 
chosen randomly if the argument has the value \"random\" (which is the 
default), or it is chosen randomly with a fixed given sum of its
entries, in  which case the argument should be a list of the form 
{\"random\",sum_value}."

PlotDualTiling::usage =
"PlotDualTiling[tiling,showlines,showpoints,linewidth,ptsize] produces
a plot of a tiling obtained from the procedure DualizeGrid; it is 
assumed that the first argument is in the form that is produced by this
routine. The other arguments are optional: showlines, True by default,
determines whether lines are plotted; showpoints, False by default,
determines whether points are plotted as circles; and finally linewidth
and ptsize fix the width of lines and size of points in the plot, their
default values being 1/200 and 1/100, respectively."

DecorateGrid::usage =
"DecorateGrid[til_,decorations_]
"

GetSlice::usage =
"GetSlice[til_,lnum_,vnum_]
"

PlotSlice::usage = 
"PlotSlice[tiling,line number,vector number,showlines,showpoints,linewidth,ptsize]
shows a slice of the tiling, specified by the Ammann line and starvector. Tiles
corresponding to the intersection points of this plane are plotted. If there are
no defects, then there will be no overlapping tiles in any given slice.
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

Unprotect[PlotGridVectors,
        DualizeGrid,
        PlotDualTiling,
		DecorateGrid,
		GetSlice,
		PlotSlice,
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

Clear[PlotGridVectors,
        DualizeGrid,
        PlotDualTiling,
		DecorateGrid,
		GetSlice,
		PlotSlice,
		GenerateRotatedTile,
		GenerateRotatedTiles,
		GetTriples,
		PlotRotatedTiles]

(*********************************)
(*   Setting global variables    *)
(*********************************)

$tau = (1+Sqrt[5])/2;
$thta0 = ArcCos[1/(2$tau-1)];(* = ArcCos[1/Sqrt[5]] for QC*)
$thta1 = ArcCos[($tau-1)/2];(* = 2Pi/5 for QC*)
$r = Union[{{0,0,1}},
				Table[CoordinateTransform["Spherical"->"Cartesian",{1,$thta0,x}],
				{x,{-2$thta1,-$thta1,0,$thta1,2$thta1}}]];
$triples = Select[Tuples@{Range[6],Range[6],Range[6]},Less@@#&];

(*********************************)
(*         Begin functions       *)
(*********************************)

PlotGridVectors[] :=
Module[{vector},
		Show[Graphics3D[Table[Arrow[{{0,0,0},vector}],
						{vector,$r}]],
			AspectRatio -> Automatic]
]


DualizeGrid[kmin_,
            kmax_,
            gamma_:"random"
			] :=
Module[{rx,rdt,
		step,cut,regch,
		kmn,kmx,gam,nj,
		i,ki,gi,kmni,kmxi,
		j,kj,gj,kmnj,kmxj,
		k,kk,gk,kmnk,kmxk,tijk,
		intersec,tiles={},tilpoints,tillines={},tilesijk={}},
			
			rx          = Table[Cross[i,j],{i,$r},{j,$r}];
			rdt         = Table[Det[{i,j,k}],{i,$r},{j,$r},{k,$r}];
			nj          = Length[$r];
			step[i_]  := step[i] = Table[If[j==i,1,0],{j,nj}];
			
			kmn        = Which[And[Head[kmin]===List,
                                 Length[kmin]===nj,
                                 Union[Map[IntegerQ,kmin]]==={True}],
                             kmin,
                             IntegerQ[kmin],
                             Table[kmin,{nj}],
                             True,
                             Print["Error: invalid argument kmin"];
                             Return[]];
			kmx        = Which[And[Head[kmax]===List,
                                 Length[kmax]===nj,
                                 Union[Map[IntegerQ,kmax]]==={True}],
                             kmax,
                             IntegerQ[kmax],
                             Table[kmax,{nj}],
                             True,
                             Print["Error: invalid argument kmax"];
                             Return[]];
			gam        = Which[gamma==="random",
                         Print["note: translation vector gamma ",
                               "is chosen randomly!"];
                         Table[Random[Real],{nj}],
                         And[Head[gamma]===List,
                             Length[gamma]==2,
                             First[gamma]==="random",
                             NumberQ[Last[gamma]]],
                         Print["note: translation vector gamma is chosen"];
                         Print["      randomly (sum of components is ",
                               gamma[[2]],")"];
                         Apply[Append[#,gamma[[2]]-Apply[Plus,#]]&,
                               {Table[Random[Real,{-0.5,0.5}],{nj-1}]}],
                         And[Head[gamma]===List,
                             Length[gamma]===nj,
                             Union[Map[NumberQ,gamma]]==={True}],
                         gamma,
                         NumberQ[gamma],
                         Table[gamma,{nj}],
                         True,
                         Print["Error: intalid argument gamma"];
                         Return[]];

			cut[p_]   := MapThread[Max,{MapThread[Min,{p,kmx+1}],kmn}];

			Do[
             gi   = gam[[i]];
             kmni = kmn[[i]];
             kmxi = kmx[[i]];
             Do[
                gj   = gam[[j]];
				kmnj = kmn[[j]];
                kmxj = kmx[[j]];
				Do[
	                gk   = gam[[k]];
      	          kmnk = kmn[[k]];
      	          kmxk = kmx[[k]];
					Do[
						intersec = 1/rdt[[i,j,k]]*	
							((ki-gi)*rx[[j,k]] + (kj-gj)*rx[[k,i]] + (kk-gk)*rx[[i,j]]);
						tijk = Table[Ceiling[N[Dot[intersec,$r[[ii]]]+
							gam[[ii]]]],{ii,nj}];
						tijk[[i]] = ki;
						tijk[[j]] = kj;
						tijk[[k]] = kk;
						AppendTo[tilesijk,{i,j,k}];
						AppendTo[tiles,Map[cut,{tijk,tijk+step[i],tijk+step[j],tijk+step[k],
											tijk+step[i]+step[j],
											tijk+step[k]+step[i],
											tijk+step[j]+step[k],
											tijk+step[i]+step[j]+step[k]}]],
						{ki,kmni,kmxi},
						{kj,kmnj,kmxj},
						{kk,kmnk,kmxk}],
					{k,j+1,nj}],
                {j,i+1,nj-1}],
             {i,nj-2}];
			
			tilpoints = Union[Flatten[tiles,1]];
			tillines  = Flatten[Table[Join[tillines,{{t[[1]],t[[2]]},{t[[1]],t[[3]]},{t[[2]],t[[5]]},{t[[3]],t[[5]]},{t[[4]],t[[6]]},
								{t[[4]],t[[7]]},{t[[6]],t[[8]]},{t[[7]],t[[8]]},{t[[1]],t[[4]]},{t[[2]],t[[6]]},{t[[3]],t[[7]]},{t[[5]],t[[8]]}}]
								,{t,tiles}],1];
			tillines  = Union[Map[Sort,tillines]];
			{nj,tilpoints,tillines,tiles,tilesijk}
]


PlotDualTiling[til_,
               showlines_:True,
               showpoints_:False,
               linewidth_:1/500,
               ptsize_:1/100] :=
Module[{add,j,plotpoints,plotlines},

			add[p_]    := Apply[Plus, p*$r];
			plotlines = 
                 If[showlines,
                    Map[Line,Map[add,til[[3]],{2}]],
                    {}];
			plotpoints = 
                 If[showpoints,
                    Map[Point,Map[add,til[[2]]]],
                    {}];
			Show[Graphics3D[Join[{Thickness[linewidth]},
								plotlines,
								{PointSize[ptsize]},
								plotpoints]],
				ViewPoint->{1.3,-2.4,2.0}*1000000,
				ImageSize->{700,700}]
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


GetSlice[til_,
		lnum_,
		vnum_]:=
Module[{rows,tnum={},tilpoints,tillines={},tiles,tilesijk,nj},

			nj=til[[1]];
			rows=Position[til[[5]],vnum][[All,1]];
			Do[
				If[til[[4,i,1,vnum]]==lnum,
					tnum=Append[tnum,i];],
			{i,rows}];
			tiles = Part[til[[4]],tnum];
			tilesijk=Part[til[[5]],tnum];
			tilpoints = Union[Flatten[tiles,1]];
			tillines  = Flatten[Table[Join[tillines,{{t[[1]],t[[2]]},{t[[1]],t[[3]]},{t[[2]],t[[5]]},{t[[3]],t[[5]]},{t[[4]],t[[6]]},
													{t[[4]],t[[7]]},{t[[6]],t[[8]]},{t[[7]],t[[8]]},{t[[1]],t[[4]]},
													{t[[2]],t[[6]]},{t[[3]],t[[7]]},{t[[5]],t[[8]]}}],
											{t,tiles}],1];
			tillines  = Union[Map[Sort,tillines]];

		{nj,tilpoints,tillines,tiles,tilesijk}
]

PlotSlice[til_,
		lnum_,
		vnum_,
		showlines_:True,
		showpoints_:False,
		linewidth_:1/500,
		ptsize_:1/100]:=
Module[{},
		PlotDualTiling[GetSlice[til,lnum,vnum],showlines,showpoints,linewidth,ptsize]
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

Protect[PlotGridVectors,
        DualizeGrid,
        PlotDualTiling,
		DecorateGrid,
		GetSlice,
		PlotSlice,
		GenerateRotatedTile,
		GenerateRotatedTiles,
		GetTriples,
		PlotRotatedTiles]

(********************************************************************)
(*                                                                  *)
(*        End of package "AperiodicTilings`GridMethod`"             *)
(*                                                                  *)
(********************************************************************)

EndPackage[]
