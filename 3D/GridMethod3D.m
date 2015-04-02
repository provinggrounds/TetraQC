(* ::Package:: *)

(********************************************************************)
(*                                                                  *)
(* :Title:       GridMethod3D                                       *)
(*                                                                  *)
(* :Authors:     Chaney Lin, Uwe Grimm                              *)
(*                                                                  *)
(* :Context:     AperiodicTilings`GridMethod`                       *)
(*                                                                  *)
(* :Version:     unofficial                                         *)
(*                                                                  *)
(* :Date:        March 31, 2015                                     *)
(*                                                                  *)
(* :Summary:     unofficial update generalizes to 3D                *)
(*               Implementation of the "grid" or "dualization"      *)
(*               method of constructing quasiperiodic tilings       *)
(*               with n-fold rotational symmetry from a regular     *)
(*               n-fold grid of equally spaces parallel lines       *)
(*                                                                  *)
(*               This package is part of a collection of            *)
(*               Mathematica programs that were originally          *)
(*               developed for a HERAEUS summer school on           *)
(*               quasicrystals held at Chemnitz University          *) 
(*               of Technology in September 1997. For more          *)
(*               information, we refer to the proceedings:          *)
(*                                                                  *)
(*               J.-B. Suck, M. Schreiber, P. Haeussler,            *)
(*               Quasicrystals: An Introduction to structure,       *)
(*               physical properties, and applications,             *)
(*               Springer, Berlin (2002)                            *)
(*                                                                  *)
(* :References:  N. G. de Bruijn,                                   *)
(*               Algebraic theory of Penrose's non-periodic         *)
(*               tilings of the plane. I & II,                      *)
(*               Indagationes mathematicae (Proc. Kon. Ned. Akad.   *)
(*               Wet. Ser. A) 84 (1991) 39-52, 53-66                *)
(*                                                                  *)
(*               V. E. Korepin, F. Gaehler, and J. Rhyner,          *)
(*               Quasiperiodic tilings: A generalized grid-         *)
(*               projection method,                                 *)
(*               Acta Cryst. A 44 (1988) 667-672                    *)
(*                                                                  *)
(* :Description: This program allows the user to construct and      *)
(*               to plot grids of equally-spaced parallel lines     *)
(*               with n-fold rotational symmetry and the n-fold     *)
(*               tilings obtained by dualization of these grids.    *)
(*               In this procedure, the grids should be regular,    *)
(*               which means that at any point no more than two     *)
(*               of the lines intersect. A check of the regularity  *)
(*               is not yet available (Version 1.0), but it is      *)
(*               meant to be added in a later version. So far,      *)
(*               there will be no error message for non-regular     *)
(*               grids, but the outcome of the dualization for      *)
(*               such grids is, more or less, unpredictable.        *)
(*                                                                  *)
(*               The grid is determined by the symmetry, the        *)
(*               number of lines in each direction, and the shifts  *)
(*               of their relative positions. These shifts can      *)
(*               also be chosen randomly in which case the          *)
(*               resulting grid may be assumed to be regular.       *) 
(*                                                                  *)
(*               The grid vectors and the grid can be plotted       *)
(*               with the commands PlotGridVectors and PlotGrid,    *)
(*               respectively. The function DualizeGrid performs    *)
(*               the dualization step, two alternative routines     *)
(*               are available which should yield the same result.  *)
(*               Finally, the function PlotDualTiling plots the     *)
(*               tiling obtained from DualizeGrid.                  *)
(*                                                                  *)
(*               Periodic approximants to quasiperiodic tilings     *)
(*               can also be constructed within the framework of    *)
(*               the grid method by an appropriate change in the    *)
(*               grid directions, this has not been implemented     *)
(*               here, but may also be added at a later stage.      *)
(*                                                                  *)
(* :Changes:     June 15, 2004 [UG] Version 1.0 -> 1.01:            *)
(*               -> Added feature PlotGridDualTiling                *)
(*               February 11, 2014 [UG] Version 1.01 -> 1.02:       *)
(*               -> Adapted for Mathematica 8                       *)
(*               -> Removed use of obsolete Graphics`Arrow`         *)
(*                  package                                         *)    
(*                                                                  *)
(* :Bugs:        No bugs known so far.                              *)
(*               Please send bug reports to:                        *)
(*               u.grimm@physics.org                                *)
(*                                                                  *)
(* :Copyright:   This package may be copied and redistributed       *) 
(*               freely by anyone, but it may not be sold           *)
(*               commercially. Any changes to the original code     *)
(*               should be documented, or preferably suggested      *)
(*               to the authors of this package to make useful      *)
(*               changes or additions available to the community.   *)
(*                                                                  *)
(* :Disclaimer:  No guarantee for correctness of the program and    *)
(*               results obtained with the program is given. Use    *)
(*               at your own risk!                                  *)
(*                                                                  *)
(********************************************************************)

BeginPackage["AperiodicTilings`GridMethod`"]

(********************************************************************)
(*                                                                  *)
(*              Usage messages for exported functions               *)
(*                                                                  *)
(********************************************************************)

PlotGridVectors::usage =
"PlotGridVectors[n] plots the grid vectors for an n-fold symmetric
grid. The argument n must be integer."

PlotGrid::usage = 
"PlotGrid[n,kmin,kmax,gamma,plotrange,length] draws an n-fold grid
for integer values of n with parallel lines in each direction at 
positions kmin to kmax. Here, kmin and kmax may either be given as a 
list of integers of length n for odd values of n, and n/2 for even 
values of n, or as a single integer in which case it will be used for 
all grid directions. The optional argument gamma determines the shift,
it can be given in three different ways: either it is provided 
explicitly as a real vector of length n for odd values of n, and n/2 
for even values of n, or it is chosen randomly if the argument has the 
value \"random\" (which is the default), or it is chosen randomly with
a fixed given sum of its entries, in  which case the argument should 
be a list of the form {\"random\",sum_value}. The optional arguments 
plotrange and length determine the PlotRange of the graph, and the 
lengths of the grid lines."

DualizeGrid::usage =
"DualizeGrid[n,kmin,kmax,gamma,check,method] calculates the tilings 
that is dual to the n-fold grid, n integer, determined by the arguments 
which have the same meaning as in the function PlotGrid. Thus, the grid
consists of  parallel lines in each direction at positions kmin to kmax. 
Here, kmin and kmax may either be given as a list of integers of length 
n for odd values of n, and n/2 for even values of n, or as a single 
integer in which case it will be used for all grid directions. The 
optional argument gamma determines the shift, it can be given in three 
different ways: either it is provided explicitly as a real vector of 
length n for odd values of n, and n/2 for even values of n, or it is 
chosen randomly if the argument has the value \"random\" (which is the 
default), or it is chosen randomly with a fixed given sum of its
entries, in  which case the argument should be a list of the form 
{\"random\",sum_value}. The optinal argument check, whose default value 
is False, is added for later convenience as an option to check the
regularity of the grid (not yet available). Finally, an optional 
argument method with default value 1 allows you to choose between two 
different implementations of the program, the default choice apparently 
performs faster."

PlotDualTiling::usage =
"PlotDualTiling[tiling,showlines,showpoints,linewidth,ptsize] produces
a plot of a tiling obtained from the procedure DualizeTiling; it is 
assumed that the first argument is in the form that is produced by this
routine. The other arguments are optional: showlines, True by default,
determines whether lines are plotted; showpoints, False by default,
determines whether points are plotted as circles; and finally linewidth
and ptsize fix the width of lines and size of points in the plot, their
default values being 1/200 and 1/100, respectively."

PlotGridDualTiling::usage = 
"PlotGridDualTiling[n,kmin,kmax,gamma,plotrange,length,showlines,
showpoints,linewidth,ptsize] draws the grid given by PlotGrid and 
plots its dual tiling given by PlotDualTiling. The parameters 
n,...,length correspond to those of PlotGrid, the parameters 
showlines,...,ptsize correspond to those of PlotDualTiling. "

(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[PlotGridVectors,
          PlotGrid,
          DualizeGrid,
          PlotDualTiling,
          PlotGridDualTiling]

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
      PlotGrid,
      DualizeGrid,
      PlotDualTiling,
      PlotGridDualTiling]
  
PlotGridVectors[tau_] :=
Module[{thta0,thta1,rj},
	thta0=ArcCos[1/(2tau-1)];(*should be arccos(1/sqrt5) for QC*)
	thta1=ArcCos[(tau-1)/2];(*should be 2Pi/5 for QC*)
	rj=Union[{{0,0,1}},Table[CoordinateTransform["Spherical"->"Cartesian",{1,thta0,x}],{x,{-2thta1,-thta1,0,thta1,2thta1}}]];
	Show[Graphics3D[Table[Arrow[{{0,0,0},r}],
                       {r,rj}]],
        AspectRatio -> Automatic]
]

	

PlotGrid[n_Integer,
         kmin_,
         kmax_,
         gamma_:"random",
         plotran_:0,
         len_:100]:= 
Module[{v,vorth,g,j,nj,kj,pr=plotran,kmn,kmx,gam},
       v[j_]     = N[{Cos[2*Pi*(j-1)/n],Sin[2*Pi*(j-1)/n]}];
       vorth[j_] = N[{-Sin[2*Pi*(j-1)/n],Cos[2*Pi*(j-1)/n]}];
       nj        = If[Mod[n,2]==0,n/2,n];
       kmn       = Which[And[Head[kmin]===List,
                             Length[kmin]===nj,
                             Union[Map[IntegerQ,kmin]]==={True}],
                         kmin,
                         IntegerQ[kmin],
                         Table[kmin,{nj}],
                         True,
                         Print["Error: invalid argument kmin"];
                         Return[]];
       kmx       = Which[And[Head[kmax]===List,
                             Length[kmax]===nj,
                             Union[Map[IntegerQ,kmax]]==={True}],
                         kmax,
                         IntegerQ[kmax],
                         Table[kmax,{nj}],
                         True,
                         Print["Error: invalid argument kmax"];
                         Return[]];
       gam       = Which[gamma==="random",
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
                         Print["Error: invalid argument gamma"];
                         Return[]];
       Print["gamma = ",gam];
       g = Table[{(kj-gam[[j]])*v[j] - len*vorth[j],
                  (kj-gam[[j]])*v[j] + len*vorth[j]},
                 {j,nj},
                 {kj,kmn[[j]],kmx[[j]]}];
       g = Map[Line,Flatten[g,1]];
       If[pr === 0,
          pr = {4{Min[kmn]-1,Max[kmx]+1},
                4{Min[kmn]-1,Max[kmx]+1}}];
       Show[Graphics[g],
            AspectRatio -> Automatic,
            PlotRange   -> pr]
       ]; 




DualizeGrid[n_Integer,
            kmin_,
            kmax_,
            gamma_:"random",
            check_:False
			] :=
Module[{tau,thta0,thta1,r,rx,rdt,
			step,cut,regch,
			kmn,kmx,gam,nj,
			i,ki,gi,kmni,kmxi,
			j,kj,gj,kmnj,kmxj,
			k,kk,gk,kmnk,kmxk,tijk,
			intersec,tiles={},tij,tilpoints,tillines={}},

			tau         = (1+Sqrt[5])/2;
			thta0       = ArcCos[1/(2tau-1)];(*should be arccos(1/sqrt5) for QC*)
			thta1       = ArcCos[(tau-1)/2];(*should be 2Pi/5 for QC*)
			r           = Union[{{0,0,1}},Table[CoordinateTransform["Spherical"->"Cartesian",{1,thta0,x}],{x,{-2thta1,-thta1,0,thta1,2thta1}}]];
			rx       = Table[Cross[i,j],{i,r},{j,r}];
			rdt = Table[Det[{i,j,k}],{i,r},{j,r},{k,r}];
			nj              = Length[r];
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

			If[check,
             Print["regularity check not yet available"];
             regch[p_] := False,
             Print["note: no check on regularity of grid will be performed"];
             Print["      result is unpredictable if there are points where"];
             Print["      more than two lines intersect"];
             regch[p_] := False];

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
						If[regch[intersec],
							Print["warning: grid is not regular!"]];
						tijk = Table[Ceiling[N[Dot[intersec,r[[ii]]]+
							gam[[ii]]]],{ii,nj}];
						tijk[[i]] = ki;
						tijk[[j]] = kj;
						tijk[[k]] = kk;
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
			{nj,tilpoints,tillines,tiles}]


PlotDualTiling[til_,
               showlines_:True,
               showpoints_:False,
               linewidth_:1/500,
               ptsize_:1/100] :=
Module[{tau,n,thta0,thta1,nj,plotpoints,plotlines,vectors,add,j},
		n       = First[til];
		tau = (1+Sqrt[5])/2;
		thta0=ArcCos[1/(2tau-1)];(*should be arccos(1/sqrt5) for QC*)
		thta1=ArcCos[(tau-1)/2];(*should be 2Pi/5 for QC*)
		vectors=Union[{{0,0,1}},Table[CoordinateTransform["Spherical"->"Cartesian",{1,thta0,x}],{x,{-2thta1,-thta1,0,thta1,2thta1}}]];
		add[p_]:= Apply[Plus, p*vectors];
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
			ImageSize->{700,700}]];


PlotGridDualTiling[n_,
                   kmin_,
                   kmax_,
                   gamma_:"random",
                   plotrange_:0,
                   length_:100,
                   showlines_:True,
                   showpoints_:False,
                   linewidth_:1/200,
                   ptsize_:1/100] := 
Module[{gam,nj},
       nj        = If[Mod[n,2]==0,n/2,n];
       gam       = Which[gamma==="random",
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
                         Print["Error: invalid argument gamma"];
                         Return[]];
Show[GraphicsColumn[{PlotGrid[n,kmin,kmax,gam,plotrange,length],
                     PlotDualTiling[DualizeGrid[n,kmin,kmax,gam]]}]]
];


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
        PlotGrid,
        DualizeGrid,
        PlotDualTiling,
        PlotGridDualTiling]

(********************************************************************)
(*                                                                  *)
(*        End of package "AperiodicTilings`GridMethod`"             *)
(*                                                                  *)
(********************************************************************)

EndPackage[]

