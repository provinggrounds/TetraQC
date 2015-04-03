(* ::Package:: *)

(********************************************************************)
(*                                                                  *)
(* :Title:       GridMethod3D                                       *)
(*                                                                  *)
(* :Authors:     Chaney Lin                                         *)
(*                                                                  *)
(* :Date:        March 31, 2015                                     *)
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

PlotSlice::usage = 
"PlotSlice[tiling,line number,vector number,showlines,showpoints,linewidth,ptsize]
shows a slice of the tiling, specified by the Ammann line and starvector. Tiles
corresponding to the intersection points of this plane are plotted. If there are
no defects, then there will be no overlapping tiles in any given slice.
"

(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[PlotGridVectors,
          DualizeGrid,
          PlotDualTiling,
		  PlotSlice]

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
	  PlotSlice]
  
PlotGridVectors[] :=
Module[{tau,thta0,thta1,rj},
		tau   = (1+Sqrt[5])/2;
		thta0 = ArcCos[1/(2tau-1)];(*should be arccos(1/sqrt5) for QC*)
		thta1 = ArcCos[(tau-1)/2];(*should be 2Pi/5 for QC*)
		rj = Union[{{0,0,1}},
					Table[CoordinateTransform["Spherical"->"Cartesian",{1,thta0,x}],
						{x,{-2thta1,-thta1,0,thta1,2thta1}}]];
		Show[Graphics3D[Table[Arrow[{{0,0,0},r}],
						{r,rj}]],
			AspectRatio -> Automatic]
]


DualizeGrid[kmin_,
            kmax_,
            gamma_:"random",
			] :=
Module[{tau,thta0,thta1,r,rx,rdt,
		step,cut,regch,
		kmn,kmx,gam,nj,
		i,ki,gi,kmni,kmxi,
		j,kj,gj,kmnj,kmxj,
		k,kk,gk,kmnk,kmxk,tijk,
		intersec,tiles={},tilpoints,tillines={},tilesijk={}},

			tau         = (1+Sqrt[5])/2;
			thta0       = ArcCos[1/(2tau-1)];(*should be arccos(1/sqrt5) for QC*)
			thta1       = ArcCos[(tau-1)/2];(*should be 2Pi/5 for QC*)
			r           = Union[{{0,0,1}},Table[CoordinateTransform["Spherical"->"Cartesian",{1,thta0,x}],{x,{-2thta1,-thta1,0,thta1,2thta1}}]];
			rx          = Table[Cross[i,j],{i,r},{j,r}];
			rdt         = Table[Det[{i,j,k}],{i,r},{j,r},{k,r}];
			nj          = Length[r];
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
						tijk = Table[Ceiling[N[Dot[intersec,r[[ii]]]+
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
Module[{tau,thta0,thta1,r,
		add,j,
		plotpoints,plotlines},

			tau         = (1+Sqrt[5])/2;
			thta0       = ArcCos[1/(2tau-1)];(*should be arccos(1/sqrt5) for QC*)
			thta1       = ArcCos[(tau-1)/2];(*should be 2Pi/5 for QC*)
			r           = Union[{{0,0,1}},Table[CoordinateTransform["Spherical"->"Cartesian",{1,thta0,x}],{x,{-2thta1,-thta1,0,thta1,2thta1}}]];
			add[p_]    := Apply[Plus, p*r];
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


PlotSlice[til_,
		lnum_,
		vnum_,
		showlines_:True,
		showpoints_:False,
		linewidth_:1/500,
		ptsize_:1/100]:=
Module[{rows,tnum={},tilpoints,tillines={},
		tau, thta0,thta1,r,add,
		plotlines,plotpoints},

			rows=Position[til[[5]],vnum][[All,1]];
			Do[
				If[til[[4,i,1,vnum]]==lnum,
					tnum=Append[tnum,i];],
			{i,rows}];
			
			tilpoints = Union[Flatten[Part[til[[4]],tnum],1]];
			tillines  = Flatten[Table[Join[tillines,{{t[[1]],t[[2]]},{t[[1]],t[[3]]},{t[[2]],t[[5]]},{t[[3]],t[[5]]},{t[[4]],t[[6]]},
													{t[[4]],t[[7]]},{t[[6]],t[[8]]},{t[[7]],t[[8]]},{t[[1]],t[[4]]},
													{t[[2]],t[[6]]},{t[[3]],t[[7]]},{t[[5]],t[[8]]}}],
											{t,Part[til[[4]],tnum]}],1];
			tillines  = Union[Map[Sort,tillines]];

			tau         = (1+Sqrt[5])/2;
			thta0       = ArcCos[1/(2tau-1)];(*should be arccos(1/sqrt5) for QC*)
			thta1       = ArcCos[(tau-1)/2];(*should be 2Pi/5 for QC*)
			r           = Union[{{0,0,1}},Table[CoordinateTransform["Spherical"->"Cartesian",{1,thta0,x}],{x,{-2thta1,-thta1,0,thta1,2thta1}}]];
			add[p_]    := Apply[Plus, p*r];
			plotlines = 
                 If[showlines,
                    Map[Line,Map[add,tillines,{2}]],
                    {}];
			plotpoints = 
                 If[showpoints,
                    Map[Point,Map[add,tilpoints]],
                    {}];
			Show[Graphics3D[Join[{Thickness[linewidth]},
								plotlines,
								{PointSize[ptsize]},
								plotpoints]],
				ViewPoint->{1.3,-2.4,2.0}*1000000,
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
		PlotSlice]

(********************************************************************)
(*                                                                  *)
(*        End of package "AperiodicTilings`GridMethod`"             *)
(*                                                                  *)
(********************************************************************)

EndPackage[]
