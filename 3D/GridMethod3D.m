(* ::Package:: *)

(********************************************************************)
(* :Title:       GridMethod3D                                       *)
(* :Authors:     Chaney Lin, based heavily on work by Uwe Grimm    *)
(* :Date:        April 2015                                         *)
(* :Summary:     unofficial update generalizes to 3D                *)
(*                                                                  *)
(********************************************************************)

BeginPackage["AperiodicTilings`GridMethod3D`"]

(********************************************************************)
(*                                                                  *)
(*              Usage messages for exported functions               *)
(*                                                                  *)
(********************************************************************)

PlotGridVectors::usage =
"PlotGridVectors[] plots the star vectors that generate the icosahedral
quasicrystal."

DualizeGrid::usage =
"DualizeGrid[kmin,kmax,gamma] calculates the icosahedral quasicrystal
tiling using the dual grid method. Initial grid consists of parallel planes
in each direction at positions kmin to kmax, possibly displaced from
origin by phase shift gamma.

Here, kmin and kmax may either be given as:
(a) a list of integers of length n for odd values of n, and n/2 for
even values of n; or
(b) a single integer in which case it will be used for all grid directions.

The optional argument gamma determines the shift, it can be given in
three different ways:
(a) either explicitly as a real vector of length n for odd values of n,
and n/2 for even values of n; or
(b) chosen randomly if the argument has the value \"random\" (which
is the default); or
(c) chosen randomly with a fixed given sum of its entries, in  which case
the argument should be a list of the form {\"random\",sum_value}.

Output is in the form (nj,tilpoints,tillines,tiles,tilesijk), where:
nj: number of star vectors
tilpoints: positions of the vertices, in terms of star vector indices
tillines: pairs of vertices defining tile lines
tiles: sets of vertices defining tiles
tilesijk: generating triplet of star vectors defining tiles

These outputs are in index form. Points can be converted to coordinates using
IndicesToCoords[tilpoints]"

GetSlice::usage =
"GetSlice[tiling, line number, vector number] returns output in the same format
as DualizeGrid, but of a slice through the tiling (given as input),
specified by the inputs line number and vector number.

Line number corresponds to the index of the
Ammann plane, perpedicular to the star vector with index specified by
vector number.

Note that line number must be within the range of the kmin and kmax
that were specified when generating tiling.

Note that vector number must be <= the number of star vectors that generated
the tiling. For icosahedral, this is 6.

Tiles corresponding to the intersection points of this plane are obtained.
If there are no defects, then there will be no overlapping tiles in the slice."

PlotDualGrid::usage =
"PlotDualGrid[tiling,showlines,showpoints,linewidth,ptsize] produces
a plot of a tiling obtained from the procedure DualizeGrid or of a slice
obtained from the procedure GetSlice. It is assumed that the first argument is
in the form that is produced by these routines.

The other arguments are optional: showlines, True by default,
determines whether lines are plotted; showpoints, False by default,
determines whether points are plotted as circles; and finally linewidth
and ptsize fix the width of lines and size of points in the plot, their
default values being 1/200 and 1/100, respectively."

PlotSlice::usage = 
"PlotSlice[tiling,line number,vector number,showlines,showpoints,linewidth,ptsize]
shows a slice of the tiling, specified by the line number and vector number."

ExportData::usage=
"ExportData[tiling,fid] exports tiling in multiple file, with ID specified by fid.

Four files are generated, containing points, lines, tiles, and triplets,
where each triplet corresponds to the star vectors generating the corresponding
tiles. Their filenames are of the form
fid_points.csv
fid_lines.csv
fid_tiles.csv
fid_ijk.csv"

IndicesToCoords::usage=
"IndicesToCoords[tilpoints] returns the coordinates corresponding to the
points specified by tilpoints."


(********************************************************************)
(*                                                                  *)
(*                   Unprotect exported functions                   *)
(*                                                                  *)
(********************************************************************)

Unprotect[PlotGridVectors,
        DualizeGrid,
		GetSlice,
        PlotDualGrid,
		PlotSlice,
		ExportData,
		IndicesToCoords]

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
		GetSlice,
        PlotDualGrid,
		PlotSlice,
		ExportData,
		IndicesToCoords]

(*********************************)
(*   Setting global variables    *)
(*********************************)

$tau = (1+Sqrt[5])/2;
$thta0 = ArcCos[1/(2$tau-1)];(* = ArcCos[1/Sqrt[5]] for QC*)
$thta1 = ArcCos[($tau-1)/2];(* = 2Pi/5 for QC*)
$r = Union[{{0,0,1}},
				Table[CoordinateTransform["Spherical"->"Cartesian",{1,$thta0,x}],
				{x,{-2$thta1,-$thta1,0,$thta1,2$thta1}}]];

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


PlotDualGrid[til_,
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


PlotSlice[til_,
		lnum_,
		vnum_,
		showlines_:True,
		showpoints_:False,
		linewidth_:1/500,
		ptsize_:1/100]:=
Module[{},
		PlotDualGrid[GetSlice[til,lnum,vnum],showlines,showpoints,linewidth,ptsize]
]

IndicesToCoords[tilpoints_]:=
Module[{coords},
			add[p_] := Apply[Plus, p*$r];
			coords = Map[add,tilpoints];
		coords
]

ExportData[til_,fid_]:= 
Module[{fpoints,flines,ftiles,fijk},
		fpoints = StringJoin[fid,"_points.csv"];
		flines = StringJoin[fid,"_lines.csv"];
		ftiles = StringJoin[fid,"_tiles.csv"];
		fijk = StringJoin[fid,"_ijk.csv"];

		Export[fpoints,til[[2]],"CSV"];
		Export[flines,til[[3]],"CSV"];
		Export[ftiles,til[[4]],"CSV"];
		Export[fijk,til[[5]],"CSV"];
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
		GetSlice,
        PlotDualGrid,
		PlotSlice,
		ExportData,
		IndicesToCoords]

(********************************************************************)
(*                                                                  *)
(*        End of package "AperiodicTilings`GridMethod3D`"           *)
(*                                                                  *)
(********************************************************************)

EndPackage[]
