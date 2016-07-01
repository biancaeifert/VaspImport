(* ::Package:: *)

(* :Context: VaspImport` *)
(* :Author: Bianca Eifert and Christian Heiliger, Theoretical Solid State Physics, Institute for Theoretical Physics, Justus Liebig University Giessen, Germany *)
(* :Summary: The VaspImport package contains an import function for file types associated with VASP (www.vasp.at). *)
(* :Copyright: Bianca Eifert and Christian Heiliger, under the MIT license *)
(* :Package Version: 1.0.0 *)
(* :Mathematica Version: 10.0 *)


If[TrueQ[$VersionNumber<10],Print["Sorry, VaspImport only works with version 10+ of the Wolfram Language."];Abort[]];


BeginPackage["VaspImport`"];
VaspImport::usage="VaspImport[file] imports all crystal structures from a file. VaspImport works with file types associated with the Vienna Ab-Initio \
Simulation Package (VASP, www.vasp.at) and it can handle CONTCAR, POSCAR, OUTCAR, XDATCAR, and vasprun.xml. The file type is detected automatically. \
Please refer to the documentation provided on the VASP website for more information about the program and its file formats. The VaspImport package may \
not work with all versions of VASP and/or all varieties of the different file types.";
Begin["`Private`"];


(*return an error for the general case; all valid file types are overloaded later on: *)
VaspImport::notfound="Could not find input structure. (Probably a missing file or an unknown file type.)";
VaspImport[___]:=(Message[VaspImport::notfound];$Failed);


VaspImport::incompleteXDATCAR="Sorry, this is an old XDATCAR which does not contain lattice vectors and therefore cannot be imported. Copy the \
lattice vector part of an appropriate CONTCAR into the top of your XDATCAR file so that the beginning of the XDATCAR looks like a CONTCAR, and try again.";


(*POSCAR / CONTCAR / CHGCAR (header only) / job script*)
VaspImport[file_String/;(FileExistsQ[file]&&StringMatchQ[FileNameTake[file],"*cnt*"|"*job*"|"*CONTCAR*"|"*POSCAR*"|"*CHG*",IgnoreCase->True])]:=
Module[{data,jso,speciesQ,chemoffset,chgoffset,lattscale,lattvec,coord,conf},

data=Import[file,"Table"];

(*offset for job scripts that write POSCARs: *)
jso=Position[data,{"cat",">","POSCAR","<<","EOF"}];
jso=If[TrueQ[jso=={}],0,jso[[1,1]]];
(*offset if chemical species are given: *)
speciesQ=MatchQ[data[[jso+6]],{_String..}|_String];
chemoffset=If[speciesQ,1,0];
(*offset for missing "Selective Dynamics" line: *)
chgoffset=If[TrueQ[Length[data[[jso+chemoffset+8]]]>0]&&MatchQ[data[[jso+chemoffset+8,1]],_?NumericQ],-1,0];

(*lattice vectors and atoms: *)
lattvec=data[[jso+3;;jso+5]];
lattscale=data[[jso+2]];
If[TrueQ[Length[lattscale]==1],lattscale=lattscale[[1]]];
If[TrueQ[Sign[lattscale]==-1],lattscale=(-lattscale/Det[lattvec])^(1/3)];
lattvec=If[TrueQ[Length[lattscale]==3],N[lattscale[[#]]*lattvec[[#]]]&/@Range[3],N[lattscale*lattvec]];
conf=Flatten[Join[ConstantArray[#,data[[jso+chemoffset+6,#]]]&/@Range[Length[data[[jso+chemoffset+6]]]]]];
coord=N[data[[jso+chemoffset+chgoffset+9;;(jso+chemoffset+chgoffset+8+Total[data[[jso+chemoffset+6]]]),1;;3]]];

(*reprojection for cartesian coordinates: *)
If[TrueQ[Length[data[[jso+chemoffset+chgoffset+8]]]>0]&&StringMatchQ[ToString[data[[jso+chemoffset+chgoffset+8,1]]],"c*"|"k*",IgnoreCase->True],
If[TrueQ[Length[lattscale]==3],coord=#*lattscale&/@coord,coord=coord*lattscale];
coord=coord.Inverse[lattvec]];

{<|
"lattice"->lattvec,
"atomcoords"->coord,
"atomtypes"->conf,
"chemical"->If[speciesQ,data[[jso+6]],{}],
"name"->StringJoin[Riffle[ToString[#]&/@data[[jso+1]]," "]],
"file"->AbsoluteFileName[file]
|>}

];


(*OUTCAR*)
VaspImport[file_String/;(FileExistsQ[file]&&StringMatchQ[FileNameTake[file],"*out*",IgnoreCase->True])]:=
Module[{data,startmarker,lattvec,cartcoord,coord,types,conf,cutpos,pos,chem,label},

data=Import[file,"Table"];
(*lattice vectors: *)
startmarker=#[[1]]&/@Position[data,"VOLUME"];
lattvec=#[[1;;3]]&/@data[[startmarker[[#]]+5;;startmarker[[#]]+7]]&/@Range[Length[startmarker]];
(*atoms and types: *)
types=data[[Select[Position[data,"ions"],#[[2]]==1&][[1,1]],5;;]];
conf=Flatten[ConstantArray[#,types[[#]]]&/@Range[Length[types]]];
cartcoord=(#[[1;;3]]&/@data[[#+2;;#+1+Length[conf]]])&/@(#[[1]]&/@Position[data,"POSITION"]);
coord=cartcoord[[#]].Inverse[lattvec[[#]]]&/@Range[Length[lattvec]];
(*chemical species: *)
cutpos=Position[data,"VRHFIN"][[1,1]];
pos=Select[Position[data,"POTCAR:"][[All,1]],#<cutpos&][[;;-2]];
chem=StringSplit[data[[#,3]],"_"][[1]]&/@pos;
(*system name: *)
label=StringJoin[Riffle[data[[Position[data,"POSCAR"][[1,1]],3;;]]," "]];

<|
"lattice"->lattvec[[#]],
"atomcoords"->coord[[#]],
"atomtypes"->conf,
"chemical"->chem,
"name"->label,
"file"->AbsoluteFileName[file]
|>&/@Range[Length[coord]]

];


(*vasprun.xml*)
VaspImport[file_String/;(FileExistsQ[file]&&StringMatchQ[FileNameTake[file],"*xml*"|"*vasprun*",IgnoreCase->True])]:=
Module[{data,lattice,lattvec,atoms,coord,types,chem,conf,label},

data=Import[file];
lattice=Cases[data,XMLElement["varray",{"name"->"basis"},x_]:>x,Infinity][[2;;-2]];
lattvec=Flatten[Cases[#,XMLElement["v",{},x_]:>ToExpression[StringSplit[x]]],1]&/@lattice;
atoms=Cases[data,XMLElement["varray",{"name"->"positions"},x_]:>x,Infinity][[2;;-2]];
coord=Flatten[Cases[#,XMLElement["v",{},x_]:>ToExpression[StringSplit[x]]],1]&/@atoms;
types=Cases[data,XMLElement["atominfo",{},x_]:>x,Infinity];
types=Cases[types,XMLElement["set",{},x_]:>x,Infinity][[1]];
types=Transpose[Partition[Flatten[Cases[types,XMLElement["c",{},x_]:>x,Infinity]],2]];
conf=ToExpression[types[[2]]];
chem=DeleteDuplicates[types[[1]]];
label=Flatten[Cases[data,XMLElement["i",{"type"->"string","name"->"SYSTEM"},x_]:>x,Infinity]][[1]];

<|
"lattice"->lattvec[[#]],
"atomcoords"->coord[[#]],
"atomtypes"->conf,
"chemical"->chem,
"name"->label,
"file"->AbsoluteFileName[file]
|>&/@Range[Length[coord]]

];


(*XDATCAR*)
(*import code assumes that different offsets as in POSCAR/CONTCAR don't occur in XDATCAR*)
VaspImport[file_String/;(FileExistsQ[file]&&StringMatchQ[FileNameTake[file],"*xdat*",IgnoreCase->True])]:=
Catch[Module[{data=Import[file,"Table"],lattvec,splitpoints,coord,conf,chem,label,anchor},

(*in old versions, the lattice vectors are missing; refuse import from such files: *)
If[MatchQ[data[[2,1]],_String]||TrueQ[data[[4]]=={"CAR"}],Message[VaspImport::incompleteXDATCAR];Throw[$Failed]];

(*distinction between the other two types of XDATCAR is VERY fragile because there are so many different XDATCAR floating around... *)
If[TrueQ[Position[data,"configuration="]=={}]||TrueQ[Length[Position[data,data[[6]]]]==1],
(*semi-old XDATCAR with only one set of lattice vectors: *)
splitpoints=Sort[Flatten[{Position[data[[8;;]],"Konfig="][[All,1]],Position[data[[8;;]],"configuration="][[All,1]],Position[data[[8;;]],{}]}]]+7;
coord=N[data[[#+1;;(#+Total[data[[7]]]),1;;3]]]&/@splitpoints;
conf=Flatten[ConstantArray[#,data[[7,#]]]&/@Range[Length[data[[7]]]]];
conf=ConstantArray[conf,Length[coord]];
lattvec=If[TrueQ[Length[data[[2]]]==3],N[data[[2,#]]*data[[2+#]]]&/@Range[3],N[data[[2,1]]*data[[3;;5]]]];
lattvec=ConstantArray[lattvec,Length[coord]];
label=If[TrueQ[Length[data[[1]]]==0],"",StringJoin[Riffle[ToString[#]&/@data[[1]]," "]]];
label=ConstantArray[label,Length[coord]];
chem=ConstantArray[data[[6]],Length[coord]];,
(*current version with multiple lattice vectors: *)
anchor=Position[data,"configuration="][[All,1]];
coord=TakeWhile[data[[#+1;;]],(Length[#]==3&&And@@NumericQ/@#)&]&/@anchor;
conf=Table[Flatten[ConstantArray[#,data[[config-1,#]]]&/@Range[Length[data[[config-1]]]]],{config,anchor}];
chem=data[[anchor-2]];
label=If[TrueQ[Length[data[[#-7]]]==0],"",StringJoin[Riffle[data[[#-7]]," "]]]&/@anchor;
lattvec=Join[{data[[2;;5]]},data[[#-6;;#-3]]&/@anchor[[1;;-2]]];
Table[lattvec[[config]]=If[TrueQ[Length[lattvec[[config,1]]]==3],N[lattvec[[config,1,#]]*lattvec[[config,1+#]]]&/@Range[3],N[lattvec[[config,1,1]]*lattvec[[config,2;;4]]]],{config,1,Length[lattvec]}];
];

<|
"lattice"->lattvec[[#]],
"atomcoords"->coord[[#]],
"atomtypes"->conf[[#]],
"chemical"->chem[[#]],
"name"->label[[#]],
"file"->AbsoluteFileName[file]
|>&/@Range[Length[coord]]

]];


End[];
EndPackage[];


(*
The MIT License (MIT)

Copyright (c) 2016 Bianca Eifert and Christian Heiliger

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*)
