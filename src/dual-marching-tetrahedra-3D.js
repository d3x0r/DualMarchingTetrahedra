// The MIT License (MIT)
//
// Copyright (c) 2020 d3x0r
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

/**
 * Marching Tetrahedra in Javascript
 *
 * Based on Unique Research
 *   - Carbon atom latice in diamond crystal; hence (DCL - Diamond Crystal Lattice)
 *
 *
 * Javascript by d3x0r
 *  - Version 1 - original 'martching tetrahedra'
 *  - Version 2 - Single Cell - Diamond Lattice
 *  - Version 3 - Single Plane - Diamond Lattice
 *  - Version 4 - Full Form - Diamond Lattice
 *  - Version DMT-3D - Full Form - Diamond Lattice
 *                    - Dual MT
 *                 - supports collision detection ?
 */


/*
	// there is a long combination of possible 3-4 bits some 50(?) of
	// which only these 6 apply in reality

this is the vertex references on the right, and which 'vert N' applies.

   0 _ _ 1
   |\   /|
    \\2//     (above page)
    | | |
     \|/
      3


 this is the line index applied.  each 'bit' of valid is this number...
	. ____0____ .
	|\         /|     01  02 03 12 23 31
	\ \_1   3_/ /
	 |  \   /  |
	  \  \ /  /
	   \  .  /   (above page)
	  2|  |  |5
	    \ 4 /
	    | | |
	     \|/
	      .

	// of all combinations of bits in crossings is 6 bits per tetrahedron,
	// these are the only valid combinations.
	// (16 bits per cell)
	const validCombinations = [
		{ val:[1,1,1,0,0,0], // valid (vert 0, 0,1,2)   0 in, 0 out
		},
		{ val:[1,1,0,0,1,1], // valid 0,1,4,5  (0-3,1-2)
		},
		{ val:[1,0,1,1,1,0], // valid 0,2,3,4  (0-2,1-3)
		},
		{ val:[1,0,0,1,0,1], // valid (vert 1, 0,3,5)
		},
		{ val:[0,1,1,1,0,1], // valid 1,2,3,5 (0-1,2,3)
		},
		{ val:[0,1,0,1,1,0], // valid (vert 2, 1,3,4 )
		},
		{ val:[0,0,1,0,1,1], // valid (vert 3, 2,4,5 )
		},
	]

*/

var DualMarchingTetrahedra3 = (function() {
const debug_ = false;
	// static working buffers

	var sizes = 0;
	const pointHolder = [null];
	const pointStateHolder = [];
	const pointMergeHolder = [[]];
	const normalHolder = [[]];
	const crossHolder = [null];
	const contentHolder = [null];
	var bits = null; // single array of true/false per cube-cell indicating at least 1 cross happened
	
	// basis cube
	const geom = [
		[0,0,0],  // bottom layer
		[1,0,0],
		[0,1,0],
		[1,1,0],
		[0,0,1],  // 5 top layer
		[1,0,1],   // 6
		[0,1,1],   // 7
		[1,1,1],   // 8
	]

	// area of a triangle scalar
	// sqrt(A)/L  is a value between 0-MagicRange.
	//   1/  0.3102016197 =  3.22370979547
	// 
	const MagicRange  =  Math.sqrt( Math.sqrt(3) / 18 );
	const MagicScalar =  1 / Math.sqrt( Math.sqrt(3) / 18 );

	// these are the 6 computed lines per cell.
	// the numbers are the indexes of the point in the computed layer map (tesselation-5tets-plane.png)
	// every other cell has a different direction of the diagonals.

// low left, front up, left diag, bot front, front diag, bottom diag, front diag
	const linesOddMin =  [ [0,2],[0,4],[2,4],  [0,1],[1,2],[1,4]  ];
	const linesEvenMin = [ [0,2],[0,4],[0,6],  [0,1],[0,3],[0,5]  ];
	const linesMin = [linesEvenMin,linesOddMin];

	// this is the running center of the faces being generated - center of the cell being computed.
	const cellOrigin = [0,0,0];

	// see next comment...
	// the order of these MUST MATCH edgeToComp order
	const vertToDataOrig = [
		[ [4,6,5,0], [3,1,5,0], [ 2,3,6,0], [6,7,5,3], [0,6,5,3] ],
		[ [ 0,2,4,1], [4,7,5,1], [6,7,4,2], [3,1,7,2], [1,2,4,7] ],
	];
	
	// these is the point orders of the tetrahedra. (first triangle in comments at top)
	// these can be changed to match the original information on the wikipedia marching-tetrahedra page.
	// the following array is modified so the result is the actual offset in the computed data plane;
	// it is computed from the above, original, array.

	// index with [odd] [tet_of_cube] [0-3 vertex]
	// result is point data rectangular offset... (until modified)
	const vertToData = [	// updated base index to resolved data cloud offset
			[ [ 2,3,6,0], [3,1,5,0], [4,6,5,0], [6,7,5,3], [0,6,5,3] ],
			[ [ 0, 2,4,1], [3,1,7,2], [6,7,4,2], [4,7,5,1], [1,2,4,7] ],
	];

	// tetrahedra edges are in this order.
	// the point indexes listed here are tetrahedra points.  (point index into vertToDataOrig points to get composite point)
	const referenceTetEdgePoints = [ [0,1] ,[0,2], [0,3], [1,2], [2,3], [3,0] ];
	
	// input is face number (0-3)  output is the triplet of points that define that face.
	// given to define relations between mating faces of tets.
	const referenceTetFacePoints = [ [ 0,1,2 ], [0,2,3], [0,3,1], [1,3,2] ];

	// input is face number (0-3) output is the three edges that bound this face.
	// (unused?  Ordering has no sense... just is what it is? )
	const referenceTetFaceEdges = [ [0,3,1], [1,2,4], [2,0,5], [3,4,5] ]; 

	const tetCount = 5; // 5 tetrahedrons per cell...

	// these are short term working variables
	// reduces number of temporary arrays created.
	const pointOutputHolder = [0,0,0];
	// face normal - base direction to add normal, scales by angle 
	const fnorm = [0,0,0];
	// used to compute cross product for angle between faces
	const tmp = [0,0,0];
	// used to compute angle between faces, are two direction vectors from a point
	const a1t = [0,0,0];
	const a2t = [0,0,0];

	// these are used for the non-geometry-helper output
	const v_cb = new THREE.Vector3();
	const v_ab = new THREE.Vector3()
	const v_normTmp = new THREE.Vector3();
	// used to compute angle between faces, are two direction vectors from a point
	const v_a1t = new THREE.Vector3();
	const v_a2t = new THREE.Vector3();

	// a tetrahedra has 6 crossing values
	// the result of this is the index into that ordered list of intersections (second triangle in comments at top)

	// indexed with [invert][face][0-1 tri/quad] [0-2]
	// indexed with 'useFace'  and that is defined above in the order of validCombinations
	const facePointIndexesOriginal = [
			[
				[[0,1,2]],    // vert 0
				[[0,1,4],[0,4,5]],
				[[0,3,4],[0,4,2]],
				[[0,5,3]],    // vert 1
				[[1,2,5],[1,5,3]],
				[[1,3,4]],    // vert 2
				[[2,4,5]]     // vert 3
			],
			[
				[[1,0,2]],    // vert 0
				[[1,0,4],[4,0,5]],
				[[3,0,4],[4,0,2]],
				[[5,0,3]],    // vert 1
				[[5,2,1],[5,1,3]],
				[[3,1,4]],    // vert 2
				[[4,2,5]]     // vert 3
			]
	];

	
   // array of child tetrahedra for collapsing in an octree of DLC cells.
   // [even/odd parent cell]
   // all cells collapse with 'even' at 0,0,0  ... 0-1  become another 0.. 2-3 -> 1 ...
   //   [tetNumberInCell]
   //      results (right now)   [x,y,z offset (n*2)+   , (exclusion),  tets at offset to include,... (1, 2 or 4 )   ]
   //   the sum of x+y+z determines whether the tet ID is 'odd' or 'even'.. (x+y+z)&1 === odd.
   //
	const childTetList = [ [
		 [ // 0 - 7 combined
			  [ 0,0,1,  /* not 3,*/  0,1,2,4  ] /* all tets but 3 - odd */
						 ,  [ 0,0,0,  0  ]
						 ,  [ 1,0,1,  0  ]
						 ,  [ 0,1,1,  0  ]
		],
		   
		[   // 1 - 7 combined
			  [ 1,0,0, /* not 2, */  0,1,3,4  ] /* all tets but 2 - odd */
						  , [ 0,0,0,  1  ]
						  , [ 1,0,1,  1  ]
						  , [ 0,1,1,  1  ]
	   ],

	   [ // 2 - 7 combined
			   [ 0,1,0, /* not 1, */ 0,2,3,4  ] /* all tets but 1 - odd */
			,[ 0,1,1,  2  ]
			,[ 1,1,1,  2  ]
			,[ 0,0,0,  2  ]
		],

		[ // 3 - 7 combined
				 [ 1,1,1, /* not 0,*/  1,2,3,4  ] /* all tets but 0 - odd */
			,[ 0,1,1,  3  ]
			,[ 1,1,0,  3  ]
			,[ 0,1,1,  3  ]
		],

		[   // 4 - 12 combined
			  [1,0,0,   2 ]
			  ,[0,0,1,   3 ]
			  ,[0,1,0,   1 ]
			  ,[1,1,1,   0 ]
			,[ 0,0,0,   3, 4  ]
			,[ 0,1,1,   1, 4  ]
			,[ 1,0,1,   2, 4  ]
			,[ 1,1,0,   0, 4  ] 
		],

		],[
			[ // 0 - 7 combined
				[ 0,0,0, /* not 3, */ 0,1,2,4  ] /* all tets but 3 - even */
							,[ 1,0,0,  0  ]
							,[ 0,1,0,  0  ]
							,[ 0,0,1,  0  ]
			],
				
			[   // 1 - 7 combined
					[ 1,0,1, /* not 2, */ 0,1,3,4  ] /* all tets but 2 - even */
								,[ 1,0,0,  1  ]
								,[ 1,1,1,  1  ]
								,[ 0,0,1,  1  ]
			],

			[ // 2 - 7 combined
					[ 0,1,1, /* not 3,*/  0,1,2,4  ] /* all tets but 1 - even */
				,[ 0,1,0,  2  ]
				,[ 0,0,1,  2  ]
				,[ 1,1,1,  2  ]
			],

			[ // 3 - 7 combined
						[ 1,1,0,  /* not 0,*/  1,2,3,4  ] /* all tets but 0 - even */
				,[ 0,1,0,  3  ]
				,[ 1,0,0,  3  ]
				,[ 1,1,1,  3  ]
			],

			[   // 4 - 12 combined
					[1,0,0,   2,4 ]
					,[0,0,1,   3,4 ]
					,[0,1,0,   1,4 ]
					,[1,1,1,   0,4 ]
				,[ 0,0,0,   3  ]
				,[ 0,1,1,   1  ]
				,[ 1,0,1,   2  ]
				,[ 1,1,0,   0  ] 
			],
		]
	]



// [odd][invert][dir][longest%3/*min*/]
const tetCentroidFacet =[  
			[ // even
		      [ [ [0,3,2,1],[0,3,2,1],[0,2,3,1],[1,3,0,2],[1,3,0,2],[0,2,3,1]]  /* unchecked! */
			  , [ [0,1,2,3],[0,1,2,3],[1,0,2,3],[1,3,0,2],[1,3,0,2],[1,0,2,3]]
			  , [ [0,3,2,1],[1,0,3,2],[3,1,0,2],[0,3,1,2],[0,3,1,2],[1,0,2,3]]
			  , [ [1,2,3,0],[1,2,3,0],[2,3,1,0],[2,0,3,1],[1,2,0,3],[2,3,1,0]]
			  , [ [0,3,2,1],[0,3,2,1],[3,2,0,1],[3,0,2,1],[3,0,2,1],[3,2,0,1]]
			  , [ [2,3,0,1],[2,3,0,1],[3,2,0,1],[3,0,2,1],[3,0,2,1],[3,2,0,1]]
			  , [ [3,2,1,0],[3,2,1,0],[3,2,0,1],[2,0,3,1],[2,0,3,1],[1,3,2,0]]
			  , [ [2,3,0,1],[2,3,0,1],[3,2,0,1],[2,0,3,1],[1,2,0,3],[3,2,0,1]]
			  ],
			 [  [ [2,1,0,3],[1,0,3,2],[0,2,3,1],[0,2,1,3],[0,2,1,3],[0,2,3,1]]  /* a few missing checks */
			  , [ [1,0,2,3],[1,0,2,3],[1,0,2,3],[0,2,1,3],[0,2,1,3],[1,0,2,3]]
			  , [ [1,0,3,2],[1,0,3,2],[3,1,0,2],[1,0,3,2],[1,0,3,2],[1,0,2,3]]
			  , [ [2,3,1,0],[0,2,3,1],[2,3,1,0],[2,3,1,0],[0,2,3,1],[2,3,1,0]]
			  , [ [0,3,2,1],[0,3,2,1],[3,0,2,1],[3,0,2,1],[3,0,2,1],[3,0,2,1]]
			  , [ [3,0,2,1],[3,0,2,1],[3,0,2,1],[3,0,2,1],[3,0,2,1],[3,0,2,1]]
			  , [ [3,2,1,0],[3,2,1,0],[0,3,2,1],[3,2,1,0],[3,2,1,0],[0,3,2,1]]
			  , [ [0,2,3,1],[0,2,3,1],[0,2,3,1],[0,2,3,1],[0,2,3,1],[0,2,3,1]]
			  ],
			],
		  [ 
			[ [ [2,3,0,1],[3,0,1,2],[0,2,3,1],[2,0,3,1],[0,3,1,2],[0,2,3,1]]  /* unchecked! */
			, [ [0,3,2,1],[0,3,2,1],[0,2,3,1],[2,0,3,1],[0,3,1,2],[1,0,2,3]]
			, [ [3,0,1,2],[3,0,1,2],[0,2,3,1],[1,3,0,2],[0,2,1,3],[1,0,2,3]] 
			, [ [0,3,2,1],[0,3,2,1],[2,3,1,0],[3,0,2,1],[0,2,1,3],[2,3,1,0]]
			, [ [2,3,0,1],[2,3,0,1],[3,2,0,1],[2,0,3,1],[2,0,3,1],[0,1,3,2]]
			, [ [1,0,3,2],[1,0,3,2],[3,2,0,1],[2,0,3,1],[2,0,3,1],[0,1,3,2]]
			, [ [1,2,3,0],[1,2,3,0],[3,2,0,1],[3,0,2,1],[2,1,3,0],[3,2,0,1]]
			, [ [1,0,3,2],[1,0,3,2],[3,2,0,1],[1,3,0,2],[0,2,1,3],[2,0,1,3]]
			],
		   [  [ [2,1,0,3],[1,0,3,2],[0,2,3,1],[1,2,0,3],[1,2,0,3],[0,2,3,1]]  /* a few missing checks */
			, [ [2,1,0,3],[2,1,0,3],[1,0,2,3],[3,1,2,0],[1,2,0,3],[1,0,2,3]]
			, [ [1,0,3,2],[1,0,3,2],[1,0,2,3],[1,0,3,2],[1,0,3,2],[1,0,2,3]]
			, [ [2,1,0,3],[0,3,2,1],[2,3,1,0],[2,3,1,0],[0,2,3,1],[2,3,1,0]]
			, [ [3,0,2,1],[3,0,2,1],[3,0,2,1],[2,0,3,1],[2,0,3,1],[3,0,2,1]]
			, [ [0,3,2,1],[0,3,2,1],[3,0,2,1],[2,0,3,1],[2,0,3,1],[3,0,2,1]]
			, [ [3,2,1,0],[3,2,1,0],[0,3,2,1],[3,2,1,0],[3,2,1,0],[0,3,2,1]]
			, [ [0,3,2,1],[0,3,2,1],[0,2,3,1],[0,3,2,1],[0,2,3,1],[0,2,3,1]]
			],
		  ],
		 ];
	


	// these are bits that are going to 0.
	const tetMasks = [ [ 4|0, 1, 2, 4|3, 4|3 ], [ 0, 4|1, 4|2, 3, 4|3 ] ];

//----------------------------------------------------------
//  This is the real working fucntion; the above is just
//  static data in the function context; instance data for this function.
	return function(data,dims, opts) {

	const elements = opts.elements;
	var vertices = opts.vertices || []
	, faces = opts.faces || [];
	var smoothShade = opts.smoothShade || false;
	var newData = [];
	const showGrid = opts.showGrid;

	const dim0 = dims[0];
	const dim1 = dims[1];
	const dim2 = dims[2];

	// these are the possible edges formed from each tet to each other near tet.
	// [odd] [tet] [face]
	const referenceTetEdgeFormations =[ [  [ [0, (1*dim1*dim0)*tetCount+0], [0,(-1*dim0)*tetCount+2], [0,(-1)*tetCount+3],   [0,4] ]
					     , [ [1, ( 1 )*tetCount+0] [1,4], [1,(-1*dim0*dim1)*tetCount+3], [1,(-1*dim0)*tetCount+1] ]
					    , [ [2, (1*dim0)*tetCount+0], [2,(-1) +1], [2,(-1*dim0*dim1)*tetCount+2], [2,4] ]
					    , [ [3, (1*dim0*dim1)*tetCount+1], [3,4], [3,(1*dim0)*tetCount+3], [3,( 1 )*tetCount+2] ]
					    , [ [4, 0], [4,1], [4,2], [4,3] ]
					    ],
					    [ [ [ 0, (-1)*tetCount+1], [0, (-1*dim0)*tetCount+2], [0,(-1*dim0*dim1)*tetCount+0], [0,4] ]
					    , [ [ 1, ( 1 )*tetCount+2], [1, (1*dim0)*tetCount+1], [1, (1*dim0*dim1)*tetCount+3], [1,4] ]
					    , [ [ 2, (1*dim0*dim1)*tetCount+2], [2, (-1)*tetCount+3], [2, (1*dim0)*tetCount+0], [2,4] ]
					    , [ [ 3, (1*dim0*dim1)*tetCount+1], [3, (-1*dim0)*tetCount+3], [3, 4], [3,( 1 )*tetCount+0] ]
					    , [ [4, 0], [4,1], [4,2], [4,3] ]
					    ]
					];


//	const referenceTetEdgeFormationsAroundTetEdges=[ [ [0],[1] ], [ [0],[2] ], [ [0], [3] ], [ [0], [4] ], [ [ 1],[2] ], [ [ 1], [3] ], [ [  2], [3]] 
//						]
	// this is, for even/odd, for each tet, for each edge, the other mating tet index (plus data offset)

	const TetFaceMapping = [ [  [ ( 1 *dim0*dim1 ) * tetCount    + 0, (-1*dim0)*tetCount + 2, (-1)*tetCount           + 3, (0)*tetCount       + 4]
			            ,  [ ( 1 ) * tetCount            + 0, (0)*tetCount       + 4, (-1*dim0*dim1)*tetCount + 3, (-1*dim0)*tetCount + 1]
			            ,  [ ( 1 *dim0 ) * tetCount      + 0, (-1)*tetCount      + 1, (-1*dim0*dim1)*tetCount + 2, (0)*tetCount       + 4]
			            ,  [ ( 1 *dim0*dim1 ) * tetCount + 1, (0)*tetCount       + 4, (1*dim0)*tetCount       + 3, ( 1 )*tetCount       + 2]
			            ,  [ ( 0 ) * tetCount            + 0, (0)*tetCount       + 3, (0)*tetCount            + 2, (0)*tetCount       + 1]
			            ]
					/* Duplicated from above - this will be a bad map... */
			          , [  [ ( -1 ) * tetCount           + 1, (-1*dim0)*tetCount + 2, (-1*dim0*dim1)*tetCount + 0, (0)*tetCount       + 4]
			            ,  [ ( 1 ) * tetCount            + 2, (1*dim0)*tetCount  + 1, (-1*dim0*dim1)*tetCount + 3, (-1*dim0)*tetCount + 3]
			            ,  [ ( 1 *dim0*dim1 ) * tetCount + 2, (-1)*tetCount      + 3, (1*dim0)*tetCount       + 1, (0)*tetCount       + 4]
			            ,  [ ( 1 *dim0*dim1 ) * tetCount + 1, (-1*dim0)*tetCount + 3, (0)*tetCount            + 4, ( 1 )*tetCount       + 0]
				    ,  [ ( 0 ) * tetCount            + 0, (0)*tetCount       + 3, (0)*tetCount            + 1, (0)*tetCount       + 2]
			            ]	
			          ];
	


	// indexed with [odd]
	//   [tet] [edge] 
	const edgeToComp = [
		[ [1*dim0*dim1*6  + 0,1*dim0*dim1*6  + 3,1,1*dim0*dim1*6  + 4,5,2]
		 ,[1*6 + 0,1*6 + 2,4,1*6 + 1,5,3]
		 ,[1*dim0*6  + 3,1*dim0*6  + 1,0,1*dim0*6  + 5,2,4]
		 ,[1*dim0*dim1*6  + 1*dim0*6  + 3,1*dim0*dim1*6  + 4,1*dim0*6  + 5,1*dim0*dim1*6  + 1*6 + 0,1*6 + 2,1*dim0*6  + 1*6 + 1]
		 ,[2,5,4,1*dim0*dim1*6  + 4,1*6 + 2,1*dim0*6  + 5]]
		,[[0,1,3,2,5,4]
		 ,[1*dim0*dim1*6  + 4,1*dim0*dim1*6  + 3,5,1*dim0*dim1*6  + 1*6 + 0,1*6 + 1,1*6 + 2]
		 ,[1*dim0*dim1*6  + 1*dim0*6  + 3,1*dim0*dim1*6  + 0,1*dim0*6  + 1,1*dim0*dim1*6  + 4,2,1*dim0*6  + 5]
		 ,[1*6 + 0,1*dim0*6  + 1*6 + 1,1*dim0*6  + 3,1*6 + 2,1*dim0*6  + 5,4]
		 ,[4,5,1*6 + 2,2,1*dim0*dim1*6  + 4,1*dim0*6  + 5] // center.
		 ]
		];


	meshCloud( data,dims );
	return null;

function makeList() {
		var context_stack = {
			first : null,
			last : null,
			saved : null,
			push(node,q) {
				var recover = this.saved;
				if( recover ) { this.saved = recover.next; recover.node = node; recover.next = null; recover.prior = this.last; }
				else { recover = { node : node, q:q, next : null, prior : this.last }; }
				if( !this.last ) this.last = this.first = recover;
				else { let x, x_=null; for( x = this.first; x; x_=x, x=x.next ) { if( recover.q > x.q ) { 
					// largest q is last in stack.
					if( !x_ ) { recover.next = this.first; this.first.prior=recover; this.first = recover; break; }
					else { x_.next = recover; recover.prior = x_; recover.next = x; x.prior = recover; break;}
				} if( !x ) { if( x_ ) x_.next = this.last = recover; } } }
				this.length++;
			},
			pop() {
				var result = this.last;
				if( !result ) return null;
				if( !(this.last = result.prior ) ) this.first = null;
				result.next = this.saved; this.saved = result;
				this.length--;
				return result.node;
			},
			length : 0,

		};
	return context_stack;

}



function PointState(v,type1,type2,typeDelta) {
	let result = {
		id:pointStateHolder.length,
		normalBuffer:[0,0,0],
		vertBuffer:[v[0],v[1],v[2]],
		visits:0,
		type1:type1,
		type2:type2,
		typeDelta:typeDelta, // saved just as meta for later
        }
	pointStateHolder.push(result );
	return result;
}

// types is an array of 6 elements.  0-1,  2-3, 4-5  are the textures at the ends 
// of the edges that p1, p2, p3 are on respectivly.
// deltas is an array of 3 floats that is the scale through that map.


function TetVert( p, n, p1, p2, p3, types, deltas ) {
	var del1 = [p2.normalBuffer[0]-p1.normalBuffer[0],p2.normalBuffer[1]-p1.normalBuffer[1],p2.normalBuffer[2]-p1.normalBuffer[2] ];
	var del2 = [p3.normalBuffer[0]-p2.normalBuffer[0],p3.normalBuffer[1]-p2.normalBuffer[1],p3.normalBuffer[2]-p2.normalBuffer[2] ];
	var del3 = [p1.normalBuffer[0]-p3.normalBuffer[0],p1.normalBuffer[1]-p3.normalBuffer[1],p1.normalBuffer[2]-p3.normalBuffer[2] ];

	var l1 = del1[0]+del1[1] +del1[2];
	var l2 = del2[0]+del2[1] +del2[2];
	var l3 = del3[0]+del3[1] +del3[2];
	//var l1 = Math.sqrt( del1[0]*del1[0]+del1[1]*del1[1] +del1[2]*del1[2] );
	//var l2 = Math.sqrt( del2[0]*del2[0]+del2[1]*del2[1] +del2[2]*del2[2] );
	//var l3 = Math.sqrt( del3[0]*del3[0]+del3[1]*del3[1] +del3[2]*del3[2] );
	var l = l1+l2+l3;

	// bigger than the grey circle; which puts them all co-linear.
	if( l > 0.1 ) {

		var a = [0,0,0];
		cross(a, del1, del2 );
		const asqr = Math.sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])/2
		const nA = asqr * MagicScalar;
		if( nA < 0.1 ) { // normalized relative area 0-0.6 is a sliver of a triangle... above this is a reasonable spread.
			// there's barely even a surface here to have a point.
			// the normals are basically all in the same plane, but not parallel.
			// plane intersections would result in at most 3 parallel lines between the 3 planes.
			p[0] = (p1.vertBuffer[0]+p2.vertBuffer[0]+p3.vertBuffer[0])/3;
			p[1] = (p1.vertBuffer[1]+p2.vertBuffer[1]+p3.vertBuffer[1])/3;
			p[2] = (p1.vertBuffer[2]+p2.vertBuffer[2]+p3.vertBuffer[2])/3;
			n[0] = p1.normalBuffer[0];
			n[1] = p1.normalBuffer[1];
			n[2] = p1.normalBuffer[2];
			
		}else {
			// intersection of planes with normals specified should result in 2 point...
			// no matter which point we choose; but we choose to calculate using the
			// two longest sides for the cross product.
			const tmp = [0,0,0];
			const tmp2 = [0,0,0];
			const t = [0,0];
			var use_n1 = p1;
			var use_n2 = p2;
			var use_n3 = p3;
			if( l1 < l2 ) {  //  p2->p1 ratio vs perimeter - want the shortest opposing side.
				if( l1 < l3 ) { // p3->p2 ratio vs perimeter - want the shortest opposing side.
				        // shortest between 1,2.
					// n3 is the furthest.
					// intersect planes.
					// 1 is least - use p3->p2,  p1->p3
					use_n1 = p3;
					use_n2 = p1;
					use_n3 = p2;
				}else
				{ 
					// 3 is least - use p1->p2, p3->p2
					// 	
					use_n1 = p2;
					use_n2 = p3;
					use_n3 = p1;
				}
			} else {
				if( l2 < l3 ) { // p3->p2 ratio vs perimeter - want the shortest opposing side.
					// 2 is least - use p2->p1, p1->p3
				}else
				{ 
					// 3 is least - use p1->p2, p3->p2
					// 	
					use_n1 = p2;
					use_n2 = p3;
					use_n3 = p1;
				}
			} 
			cross( tmp, use_n1.normalBuffer, use_n2.normalBuffer );  // slope of line intersecting these
			cross( tmp2, use_n1.normalBuffer, tmp ); // normal of the plane of N3-N2
			IntersectLineWithPlane( tmp2, use_n2.vertBuffer, use_n1.normalBuffer, use_n1.vertBuffer, t );
			tmp2[0] = use_n1.vertBuffer[0]+use_n1.normalBuffer[0]*t[1];
			tmp2[1] = use_n1.vertBuffer[1]+use_n1.normalBuffer[1]*t[1];
			tmp2[3] = use_n1.vertBuffer[2]+use_n1.normalBuffer[2]*t[1];
			IntersectLineWithPlane( tmp, tmp2, use_n3.normalBuffer, use_n3.vertBuffer, t );
			p[0] = use_n3.vertBuffer[0]+use_n3.normalBuffer[0]*t[1];
			p[1] = use_n3.vertBuffer[1]+use_n3.normalBuffer[1]*t[1];
			p[2] = use_n3.vertBuffer[2]+use_n3.normalBuffer[2]*t[1];

			// scale each normal, by the inner angle of this face's vert...
			ScaleNormalsByAngle( n, use_n1.normalBuffer, use_n1.vertBuffer, use_n2.vertBuffer, use_n3.vertBuffer );
			ScaleNormalsByAngle( n, use_n2.normalBuffer, use_n2.vertBuffer, use_n3.vertBuffer, use_n1.vertBuffer );
			ScaleNormalsByAngle( n, use_n3.normalBuffer, use_n3.vertBuffer, use_n1.vertBuffer, use_n2.vertBuffer );

			// linear naive average  normals
			//n[0] = use_n1.normalBuffer[0]+use_n2.normalBuffer[0]+use_n3.normalBuffer[0]
			//n[1] = use_n1.normalBuffer[1]+use_n2.normalBuffer[1]+use_n3.normalBuffer[1]
			//n[2] = use_n1.normalBuffer[2]+use_n2.normalBuffer[2]+use_n3.normalBuffer[2]
			const nLen = 1/Math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
			n[0] *= nLen;
			n[1] *= nLen;
			n[2] *= nLen;
			
		}
	}else {
		// all normals are co-linear.... 
		// all should also be in the same direction, and normal to the face; it would be hard to imagine computing 
		// a surface with a skewed/sheared normal.
		p[0] = (p1.vertBuffer[0]+p2.vertBuffer[0]+p3.vertBuffer[0])/3;
		p[1] = (p1.vertBuffer[1]+p2.vertBuffer[1]+p3.vertBuffer[1])/3;
		p[2] = (p1.vertBuffer[2]+p2.vertBuffer[2]+p3.vertBuffer[2])/3;
		n[0] = p1.normalBuffer[0];
		n[1] = p1.normalBuffer[1];
		n[2] = p1.normalBuffer[2];
	}

}




function IntersectLineWithPlane( Slope, Origin,  // line m, b
					 n, o,  // plane n, o
					 r ) {
		var a,b,c,cosPhi, t; // time of intersection
		a = ( Slope.x * n.x + Slope.y * n.y + Slope.z * n.z );
		r.x = 0.0; r.y = 0.0;
	        
		if( a == 0.0 ) return;
	        
		b = length( Slope );
		c = length( n );
		if( b == 0.0 || c == 0.0 ) return; // bad vector choice - if near zero length...
	        
		cosPhi = a / ( b * c );
		t = ( n.x * ( o.x - Origin.x ) + n.y * ( o.y - Origin.y ) + n.z * ( o.z - Origin.z ) ) / a;
	        
		if( cosPhi > 0.0 || cosPhi < 0.0 ) { // at least some degree of insident angle
			r.x = cosPhi; r.y = t;
		} else
			return;
	}

function cross(o,a,b) {
	o[0] = a[1]*b[2]-a[2]*b[1];
	o[1] = a[2]*b[0]-a[0]*b[2];
	o[2] = a[0]*b[1]-a[1]*b[0];
	return o;
}

function ncross(o,a,b) {
	const an = [0,0,0];
	const bn = [0,0,0];
	const al = 1/Math.sqrt(a[1]*a[1]+a[0]*a[0]+a[2]*a[2]);
	const bl = 1/Math.sqrt(b[1]*b[1]+b[0]*b[0]+b[2]*b[2]);
	an[0] = a[0]*al;
	an[1] = a[1]*al;
	an[2] = a[2]*al;

	bn[0] = b[0]*bl;
	bn[1] = b[1]*bl;
	bn[2] = b[2]*bl;

	o[0] = an[1]*bn[2]-an[2]*bn[1];
	o[1] = an[2]*bn[0]-an[0]*bn[2];
	o[2] = an[0]*bn[1]-an[1]*bn[0];
	return o;
}

function IntersectNormals( po, no, p1, n1, p2, n2, p3, n3 ) {

	const tmp = [0,0,0];
	const tmp2 = [0,0,0];
	const t = [0,0];

	if( (n1[0]-n2[0]) + (n1[1]-n2[1]) + (n1[2]-n2[2]) < 0.1 ) {
		// parallel.
		if( (n1[0]-n3[0]) + (n1[1]-n3[1]) + (n1[2]-n3[2]) < 0.1 ) {
			// all 3 are parallel... return the average point.
			po[0] = (p1[0]+	p2[0]+ p3[0])/3;
			po[1] = (p1[1]+	p2[1]+ p3[1])/3;
			po[2] = (p1[2]+	p2[2]+ p3[2])/3;
			no[0] = n1[0];
			no[1] = n1[1];
			no[2] = n1[2];
		} else {
			// flex between N1 and N3
			// vertical plane between N1 and N2
		
			cross( tmp, n1, n3 );
			cross( tmp2, n1, tmp );
			IntersectLineWithPlane( tmp2, p3, n1, p1, t );
			tmp2[0] = p1[0]+n1[0]*t[0];
			tmp2[1] = p1[1]+n1[1]*t[0];
			tmp2[3] = p1[2]+n1[2]*t[0];
			IntersectLineWithPlane( tmp, tmp2, n3, p3, t );
			po[0] = p3[0]+n3[0]*t[0];
			po[1] = p3[1]+n3[1]*t[0];
			po[2] = p3[2]+n3[2]*t[0];
		}

	} else {
		if( (n1[0]-n3[0]) + (n1[1]-n3[1]) + (n1[2]-n3[2]) < 0.1 ) {
			// flex between 1 and 2
			// vertical plane between 1 and 3
			cross( tmp, n1, n2 );
			cross( tmp2, n1, tmp );
			IntersectLineWithPlane( tmp2, p2, n1, p1, t );
			tmp2[0] = p1[0]+n1[0]*t[0];
			tmp2[1] = p1[1]+n1[1]*t[0];
			tmp2[3] = p1[2]+n1[2]*t[0];
			IntersectLineWithPlane( tmp, tmp2, n2, p2, t );
			po[0] = p3[0]+n3[0]*t[0];
			po[1] = p3[1]+n3[1]*t[0];
			po[2] = p3[2]+n3[2]*t[0];
		} else {
			// intersect planes.
			cross( tmp, n1, n2 );
			cross( tmp2, n1, tmp );
			IntersectLineWithPlane( tmp2, p2, n1, p1, t );
			tmp2[0] = p1[0]+n1[0]*t[0];
			tmp2[1] = p1[1]+n1[1]*t[0];
			tmp2[3] = p1[2]+n1[2]*t[0];
			IntersectLineWithPlane( tmp, tmp2, n3, p3, t );
			po[0] = p3[0]+n3[0]*t[0];
			po[1] = p3[1]+n3[1]*t[0];
			po[2] = p3[2]+n3[2]*t[0];
		}
	}
}


function ScaleNormalsByAngle( n, fnorm, v1, v2, v3 ) {
	
	a1t[0]=vB[0]-vA[0];a1t[1]=vB[1]-vA[1];a1t[2]=vB[2]-vA[2];
	a2t[0]=vC[0]-vA[0];a2t[1]=vC[1]-vA[1];a2t[2]=vC[2]-vA[2];

	let angle = 0;
	if( (a1t[0]*a1t[0]+a1t[1]*a1t[1]+a1t[2]*a1t[2] ) >0.00000001 && 
	    (a2t[0]*a2t[0]+a2t[1]*a2t[1]+a2t[2]*a2t[2] ) >0.00000001 )
		angle = 2*Math.acos( clamp((a1t[0]*a2t[0]+a1t[1]*a2t[1]+a1t[2]*a2t[2])/(Math.sqrt(a1t[0]*a1t[0]+a1t[1]*a1t[1]+a1t[2]*a1t[2])*Math.sqrt(a2t[0]*a2t[0]+a2t[1]*a2t[1]+a2t[2]*a2t[2] ) ), 1.0 ));
	n += fnorm[0]*angle;
	n += fnorm[1]*angle;
	n += fnorm[2]*angle;
	
}



function meshCloud(data, dims) {

	// values input to this are in 2 planes for lower and upper values
	const dataOffset = [ 0, 1, dim0, 1+dim0, 0 + dim0*dim1,1 + dim0*dim1,dim0 + dim0*dim1, 1+dim0 + dim0*dim1] ;

	const pointSchedule = makeList();

	// vertex paths 0-1 0-2, 0-3  1-2 2-3 3-1
	// this is the offset from dataOffset to the related value.

	// index with [odd] [tet_of_cube] [0-5 line index]
	// result is composite point data offset.
	
	for( let a = 0; a < 2; a++ ) for( let b = 0; b < 5; b++ ) for( let c = 0; c < 4; c++ ) vertToData[a][b][c] = dataOffset[vertToDataOrig[a][b][c]];


	// this is a computed lookup from facePointIndexes ([invert][output_face_type][0-1 triangle count][0-3 triangle point index]
	// it is actually edgeToComp[odd][tet][  FPI[invert][face_type][0-1][point indexes] ]
	// index with [odd][tet][invert][output_face_type][0-1 triangle count][0-3 triangle point index]
	const facePointIndexes = [ 
	];

	for( let odd=0; odd < 2; odd++) {
		let t; 
		facePointIndexes.push( t = [] )
		for( let tet = 0; tet < 5; tet++ ){
			t.push(  [
				[ [ [edgeToComp[odd][tet][facePointIndexesOriginal[0][0][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[0][0][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[0][0][0][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[0][1][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[0][1][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[0][1][0][2]]]
				, [  edgeToComp[odd][tet][facePointIndexesOriginal[0][1][1][0]],edgeToComp[odd][tet][facePointIndexesOriginal[0][1][1][1]],edgeToComp[odd][tet][facePointIndexesOriginal[0][1][1][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[0][2][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[0][2][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[0][2][0][2]]]
				, [  edgeToComp[odd][tet][facePointIndexesOriginal[0][2][1][0]],edgeToComp[odd][tet][facePointIndexesOriginal[0][2][1][1]],edgeToComp[odd][tet][facePointIndexesOriginal[0][2][1][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[0][3][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[0][3][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[0][3][0][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[0][4][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[0][4][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[0][4][0][2]]]
				, [  edgeToComp[odd][tet][facePointIndexesOriginal[0][4][1][0]],edgeToComp[odd][tet][facePointIndexesOriginal[0][4][1][1]],edgeToComp[odd][tet][facePointIndexesOriginal[0][4][1][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[0][5][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[0][5][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[0][5][0][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[0][6][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[0][6][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[0][6][0][2]]] ]
			]
			, [ [   [edgeToComp[odd][tet][facePointIndexesOriginal[1][0][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[1][0][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[1][0][0][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[1][1][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[1][1][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[1][1][0][2]]]
				, [  edgeToComp[odd][tet][facePointIndexesOriginal[1][1][1][0]],edgeToComp[odd][tet][facePointIndexesOriginal[1][1][1][1]],edgeToComp[odd][tet][facePointIndexesOriginal[1][1][1][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[1][2][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[1][2][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[1][2][0][2]]]
				, [  edgeToComp[odd][tet][facePointIndexesOriginal[1][2][1][0]],edgeToComp[odd][tet][facePointIndexesOriginal[1][2][1][1]],edgeToComp[odd][tet][facePointIndexesOriginal[1][2][1][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[1][3][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[1][3][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[1][3][0][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[1][4][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[1][4][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[1][4][0][2]]]
				, [  edgeToComp[odd][tet][facePointIndexesOriginal[1][4][1][0]],edgeToComp[odd][tet][facePointIndexesOriginal[1][4][1][1]],edgeToComp[odd][tet][facePointIndexesOriginal[1][4][1][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[1][5][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[1][5][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[1][5][0][2]]] ]
				, [ [edgeToComp[odd][tet][facePointIndexesOriginal[1][6][0][0]],edgeToComp[odd][tet][facePointIndexesOriginal[1][6][0][1]],edgeToComp[odd][tet][facePointIndexesOriginal[1][6][0][2]]] ]
			] ] )	
		}
	}


	if( dim0*dim1*dim2*6 > sizes ) {
		sizes = dim0 * dim1 * dim2*6;
		bits = new Uint8Array(dim0*dim1*dim2);
		pointHolder[0] = new Uint32Array(sizes);
		crossHolder[0] = new Uint8Array(sizes);
		contentHolder[0] = new Uint8Array(dim0*dim1*dim2*5);  // 1 bit for each tet in each cell. (5)
		for( let zz = normalHolder[0].length; zz < sizes; zz++ ) {
			normalHolder[0].push( null );
		}
		for( let zz = pointMergeHolder[0].length; zz < (dim0*dim1*dim2); zz++ ) {
			pointMergeHolder[0].push( null ); 
		}
	}
	
	// all work space has been allocated by this point.
	// now, for each layer, compute inside-outside crossings ('cross').. which interestingly relates to 'cross product' but in a 2D way...

	const points  = pointHolder[0];
	const pointMerge  = pointMergeHolder[0];
	const normals = normalHolder[0];
	const crosses = crossHolder[0];
	const content = contentHolder[0];
	for( let zero = 0; zero < dim0*dim1*dim2; zero++ ) {
		pointMerge[zero] = null;
		// make sure to reset this... 
		for( let crz = 0; crz < 6; crz++ ) crosses[zero*6+crz] = 0;
		for( let crz = 0; crz < 5; crz++ )  { content[zero*5+crz] = 0; normals[zero*5+crz] = null; }
	}
	for( var z = 0; z < dim2; z++ ) {
	
		let odd = 0;
		let zOdd = z & 1;
		cellOrigin[2] = z-0.5;

		// compute one layer (x by y) intersections (cross from inside to outside).
		// each cell individually has 16 intersections
		// the first cell needs 9 intersections computed; the subsequent cells providing the computation for the 7 'missing'
		// 3 intersections per cell after the first layer can be copied; but shift in position (moving from the top to the bottom)
		// 
		for( var y = 0; y < dim1-1; y++ ) {
			cellOrigin[1] = y-0.5;
			for( var x = 0; x < dim0-1; x++ ) {
				odd = (( x + y ) &1) ^ zOdd;
				//if( x > 5 || y > 5 || z > 5 ) continue;
				cellOrigin[0] = x-0.5;
	
				const baseHere = (x+0 + y*dim0 + z*(dim0*dim1))*6;
				const baseOffset = x+0 + y*dim0 + z * dim0*dim1;
				const lineArray = linesMin[odd];
				bits[baseOffset] = 0;
				//console.log( "Set bits to 0", baseOffset)

				for( let l = 0; l < 6; l++ ) {
					const p0 = lineArray[l][0];
					const p1 = lineArray[l][1];
	
					if( (x == (dim0-1)) &&( (p0 & 1) || (p1 &1) )) {
						// this is an overflow, need to push fake data....
						points[baseHere+l] = null;
						crosses[baseHere+l] = 0;
						continue;
					}
					if( (y == (dim1-1)) &&( (p0 & 2) || (p1 &2) )) {
						// this is an overflow, need to push fake data....
						points[baseHere+l] = null;
						crosses[baseHere+l] = 0;
						continue;
					}
					if( (z == (dim2-1)) &&( (p0 & 4) || (p1 & 4) )) {
						// this is an overflow, need to push fake data....
						points[baseHere+l] = null;
						crosses[baseHere+l] = 0;
						continue;
					}

					const data0=baseOffset+dataOffset[p0];
					const data1=baseOffset+dataOffset[p1];
	
					d=-data[data0]; e=-data[data1];
	
					if( ( d <= 0 && e >0  )|| (d > 0 && e <= 0 ) ){
						let t;
						let normal = null;
						//console.log( "x, y is a cross:", x+y*dim0,(x+y*dim0)*6, crosses.length, baseOffset+l, x, y, p0, p1, data0, data1, `d:${d} e:${e}` );
						if( e <= 0 ) {
							(t = -e/(d-e));
// --V-V-V-V-V-V-V-V-- CREATE OUTPUT POINT(VERTEX) HERE --V-V-V-V-V-V-V-V--
							normal = null;
							if( t < 0.00001 )  {
								normal = pointMerge[baseOffset+dataOffset[p1]];
							} else if( t > 1-0.001 ){
								normal = pointMerge[baseOffset+dataOffset[p0]];
							}
							if( !normal ) {
								pointOutputHolder[0] = x + geom[p1][0]+( geom[p0][0]- geom[p1][0])* t;
								pointOutputHolder[1] = y + geom[p1][1]+( geom[p0][1]- geom[p1][1])* t;
								pointOutputHolder[2] = z + geom[p1][2]+( geom[p0][2]- geom[p1][2])* t;
								normal = PointState( pointOutputHolder
									, elements.data[data0]
									, elements.data[data1]
									, t
								);
								if( t < 0.00001 )  {
									pointMerge[baseOffset+dataOffset[p1]] = normal;
								} else if( t > 1-0.001 ) {
									pointMerge[baseOffset+dataOffset[p0]] = normal;	
								}
							}
							points[baseHere+l] = normal.id;
// --^-^-^-^-^-^-- END OUTPUT POINT(VERTEX) HERE --^-^-^-^-^-^--
						} else {
							(t = -d/(e-d));
// --V-V-V-V-V-V-V-V-- OUTPUT POINT 2 HERE --V-V-V-V-V-V-V-V--
							normal = null;
							if( t < 0.0001 )  {
								normal = pointMerge[baseOffset+dataOffset[p0]];
							} else if( t > 1-0.001 )  {
								normal = pointMerge[baseOffset+dataOffset[p1]];
							}
							if( !normal ) {
								pointOutputHolder[0] = x + geom[p0][0]+( geom[p1][0]- geom[p0][0])* t;
								pointOutputHolder[1] = y + geom[p0][1]+( geom[p1][1]- geom[p0][1])* t;
								pointOutputHolder[2] = z + geom[p0][2]+( geom[p1][2]- geom[p0][2])* t;
								normal = PointState( pointOutputHolder
									, elements.data[data0]
									, elements.data[data1]
									, t
								);
								if( t < 0.0001 ) 
									pointMerge[baseOffset+dataOffset[p0]] = normal;
								else if( t > 1-0.001 )
									pointMerge[baseOffset+dataOffset[p1]] = normal;
							}
							points[baseHere+l] = normal.id;
// --^-^-^-^-^-^-- END  OUTPUT POINT 2 HERE --^-^-^-^-^-^--
						}
						//console.log( "Set Used:", x, y, z, odd, baseOffset, baseOffset*5, bits[baseOffset], l);
						crosses[baseHere+l] = 1;
						//bits[baseOffset] = 1; // set any 1 bit is set here.
					}
					else {
						//console.log( "x,y does not cross", x+y*dim0,(x+y*dim0)*6, crosses.length, baseOffset+l, x, y, p0, p1, data0, data1, `d:${d} e:${e}` ); 
						crosses[baseHere+l] = 0;
					}
				}
			}
		}
	}


	for( var z = 0; z < dim2; z++ ) {
	
		let odd = 0;
		let zOdd = z & 1;

		// for all bounday crossed points, generate the faces from the intersection points.
		for( var y = 0; y < dim1; y++ ) {
			for( var x = 0; x < dim0; x++ ) {
				let tetSkip = 0;
				const baseOffset = (x + (y*dim0) + z*dim0*dim1)*6;

				if( x >= (dim0-1)) tetSkip |= 1;
				if( y >= (dim1-1)) tetSkip |= 2;
				if( z >= (dim2-1)) tetSkip |= 4;
				const dataOffset = x + (y*dim0) + z*dim1*dim0;
	        		odd = (( x + y ) &1) ^ zOdd;
				for( tet = 0; tet < 5; tet++ ) {
					if( tetMasks[odd][tet] & tetSkip ) continue;

					let f;
					let invert = 0;
					let useFace = 0;

					// this is 'valid combinations' check.
					if( crosses[ baseOffset+edgeToComp[odd][tet][0] ] ) {
						//console.log( `Output: odd:${odd} tet:${tet} x:${x} y:${y} a:${JSON.stringify(a)}` );
						if( crosses[ baseOffset+edgeToComp[odd][tet][1] ] ) {
							if( crosses[ baseOffset+edgeToComp[odd][tet][2] ] ) {
								// lower left tet. // source point is 0
								useFace = 1;
								invert = ( data[dataOffset+vertToData[odd][tet][0]] >= 0 )?1:0;
							} else {
								if( crosses[ baseOffset+edgeToComp[odd][tet][4] ] && crosses[ baseOffset+edgeToComp[odd][tet][5] ]) {
									// source point is 2? 1?   (0?3?)
									useFace = 2;
									invert = ( data[dataOffset+vertToData[odd][tet][0]] >= 0 )?1:0 ;
								}
							}
						} else {
							if( crosses[ baseOffset+edgeToComp[odd][tet][2]] && crosses[ baseOffset+edgeToComp[odd][tet][3]] && crosses[ baseOffset+edgeToComp[odd][tet][4] ] ) {
								// source point is ? 1? 3?   (0? 2?)
								useFace = 3;
								invert = ( data[dataOffset+vertToData[odd][tet][0]] >= 0 )?1:0  ;
							}else if( crosses[ baseOffset+edgeToComp[odd][tet][3]] && crosses[ baseOffset+edgeToComp[odd][tet][5] ] ) {
								// source point is 1
								useFace = 4;
								invert = ( data[dataOffset+vertToData[odd][tet][1]] >= 0 )?1:0
							}
						}
					} else {
						if( crosses[ baseOffset+edgeToComp[odd][tet][1] ] ) {
							if( crosses[ baseOffset+edgeToComp[odd][tet][2] ] && crosses[ baseOffset+edgeToComp[odd][tet][3] ] && crosses[ baseOffset+edgeToComp[odd][tet][5] ]) {
								// 0?1?   2?3?
								useFace = 5;
								invert = ( data[dataOffset+vertToData[odd][tet][0]] >= 0 )  ?1:0
							} else if( crosses[ baseOffset+edgeToComp[odd][tet][3]] && crosses[ baseOffset+edgeToComp[odd][tet][4] ] ) {
								// source point is 2
								useFace = 6;
								invert = ( data[dataOffset+vertToData[odd][tet][2]] >= 0 ) ?1:0
							}
						} else {
							if( crosses[ baseOffset+edgeToComp[odd][tet][2] ] && crosses[ baseOffset+edgeToComp[odd][tet][4]] && crosses[ baseOffset+edgeToComp[odd][tet][5] ] ) {
								// source point is 3
								useFace = 7;
								invert = ( data[dataOffset+vertToData[odd][tet][3]] >= 0 ) ?1:0
							} else {
							}
						}
					}

					//if( useFace > 5 || useFace < 5 ) continue;
					if( useFace-- ) {
						const fpi = facePointIndexes[odd][tet][invert][useFace];
						for( var tri=0;tri< fpi.length; tri++ ){
							const ai = points[baseOffset+fpi[tri][0]];
							const bi = points[baseOffset+fpi[tri][1]];
							const ci = points[baseOffset+fpi[tri][2]];
							// ai, bi, ci are indexes into computed pointcloud layer.
// --V-V-V-V-V-V-V-V-- GENERATE POINT NORMALS --V-V-V-V-V-V-V-V--

								//  https://stackoverflow.com/questions/45477806/general-method-for-calculating-smooth-vertex-normals-with-100-smoothness
								// suggests using the angle as a scalar of the normal.
								
								// a - b - c    c->b a->b
								const vA = pointStateHolder[ai].vertBuffer;
								const vB = pointStateHolder[bi].vertBuffer;
								const vC = pointStateHolder[ci].vertBuffer;
								let v1, v2, v3;
								const AisB =  ( ( vA[0] === vB[0] ) && ( vA[1] === vB[1]  ) && ( vA[2] === vB[2]  ) );
								const AisC =  ( ( vA[0] === vC[0] ) && ( vA[1] === vC[1]  ) && ( vA[2] === vC[2]  ) );
								const BisC =  ( ( vB[0] === vC[0] ) && ( vB[1] === vC[1]  ) && ( vB[2] === vC[2]  ) );
								if( AisB || BisC || AisC ) {
								   //console.log( "zero size tri-face")
								   continue;
								}
								{
									v1 = vC;
									v2 = vB;
									v3 = vA;
								}
								if( AisC ) {
									v1 = vB;
									v2 = vC;
									v3 = vA;
								}
								if( BisC ) {
									v1 = vA;
									v2 = vC;
									v3 = vB;
								}

								//if( !vA || !vB || !vC ) debugger;
								fnorm[0] = v2[0]-v1[0];fnorm[1] = v2[1]-v1[1];fnorm[2] = v2[2]-v1[2];
								tmp[0] = v3[0]-v1[0];tmp[1] = v3[1]-v1[1];tmp[2] = v3[2]-v1[2];
								let a=fnorm[0], b=fnorm[1];
								fnorm[0]=fnorm[1]*tmp[2] - fnorm[2]*tmp[1];
								fnorm[1]=fnorm[2]*tmp[0] - a       *tmp[2];
								fnorm[2]=a       *tmp[1] - b       *tmp[0];
								let ds;
								if( (ds=fnorm[0]*fnorm[0]+fnorm[1]*fnorm[1]+fnorm[2]*fnorm[2]) > 0.00000001 ){
									ds = 1/Math.sqrt(ds);
									fnorm[0] *= ds;fnorm[1] *= ds;fnorm[2] *= ds;
								}else {
									//console.log( "1Still not happy...", fnorm, ds,vA, vB, vC );
									// b->A  c->A
									fnorm[0] = vB[0]-vA[0];fnorm[1] = vB[1]-vA[1];fnorm[2] = vB[2]-vA[2];
									tmp[0] = vC[0]-vA[0];tmp[1] = vC[1]-vA[1];tmp[2] = vC[2]-vA[2];
									let a=fnorm[0];
									fnorm[0]=fnorm[1]*tmp[2] - fnorm[2]*tmp[1];
									fnorm[1]=fnorm[2]*tmp[0] - a       *tmp[2];
									fnorm[2]=a       *tmp[1] - b       *tmp[0];
									let ds2;
									if( (ds2=fnorm[0]*fnorm[0]+fnorm[1]*fnorm[1]+fnorm[2]*fnorm[2]) > 0.00000001 ){
										ds2 = 1/Math.sqrt(ds2);
										fnorm[0] *= ds2;fnorm[1] *= ds2;fnorm[2] *= ds2;
									} else {
										//console.log( "2Still not happy...", ds2, vA, vB, vC );
										// B->C  A->C
										fnorm[0] = vA[0]-vC[0];fnorm[1] = vA[1]-vC[1];fnorm[2] = vA[2]-vC[2];
										tmp[0] = vB[0]-vC[0];tmp[1] = vB[1]-vC[1];tmp[2] = vB[2]-vC[2];
										let a=fnorm[0];
										fnorm[0]=fnorm[1]*tmp[2] - fnorm[2]*tmp[1];
										fnorm[1]=fnorm[2]*tmp[0] - a       *tmp[2];
										fnorm[2]=a       *tmp[1] - b       *tmp[0];
										let ds3;
										if( (ds3=fnorm[0]*fnorm[0]+fnorm[1]*fnorm[1]+fnorm[2]*fnorm[2]) > 0.00000001 ){
											ds3 = 1/Math.sqrt(ds3);
											fnorm[0] *= ds3;fnorm[1] *= ds3;fnorm[2] *= ds3;
										} 
										//else 
										//	console.log( "3Still not happy...", ds, vA, vB, vC );
									}
								}
	
								{
									a1t[0]=vB[0]-vA[0];a1t[1]=vB[1]-vA[1];a1t[2]=vB[2]-vA[2];
									a2t[0]=vC[0]-vA[0];a2t[1]=vC[1]-vA[1];a2t[2]=vC[2]-vA[2];

									let angle = 0;
									if( (a1t[0]*a1t[0]+a1t[1]*a1t[1]+a1t[2]*a1t[2] ) >0.00000001 && 
									    (a2t[0]*a2t[0]+a2t[1]*a2t[1]+a2t[2]*a2t[2] ) >0.00000001 )
										angle = 2*Math.acos( clamp((a1t[0]*a2t[0]+a1t[1]*a2t[1]+a1t[2]*a2t[2])/(Math.sqrt(a1t[0]*a1t[0]+a1t[1]*a1t[1]+a1t[2]*a1t[2])*Math.sqrt(a2t[0]*a2t[0]+a2t[1]*a2t[1]+a2t[2]*a2t[2] ) ), 1.0 ));
									pointStateHolder[ai].normalBuffer[0] += fnorm[0]*angle;
									pointStateHolder[ai].normalBuffer[1] += fnorm[1]*angle;
									pointStateHolder[ai].normalBuffer[2] += fnorm[2]*angle;
								}

								{
									a1t[0]=vC[0]-vB[0];a1t[1]=vC[1]-vB[1];a1t[2]=vC[2]-vB[2];
									a2t[0]=vA[0]-vB[0];a2t[1]=vA[1]-vB[1];a2t[2]=vA[2]-vB[2];
									let angle = 0;
									if( (a1t[0]*a1t[0]+a1t[1]*a1t[1]+a1t[2]*a1t[2] ) >0.00000001 && 
									    (a2t[0]*a2t[0]+a2t[1]*a2t[1]+a2t[2]*a2t[2] ) >0.00000001 ) {
											angle = 2*Math.acos( clamp((a1t[0]*a2t[0]+a1t[1]*a2t[1]+a1t[2]*a2t[2])/(Math.sqrt(a1t[0]*a1t[0]+a1t[1]*a1t[1]+a1t[2]*a1t[2])*Math.sqrt(a2t[0]*a2t[0]+a2t[1]*a2t[1]+a2t[2]*a2t[2] ) ), 1.0) );
									}

									pointStateHolder[bi].normalBuffer[0] += fnorm[0]*angle;
									pointStateHolder[bi].normalBuffer[1] += fnorm[1]*angle;
									pointStateHolder[bi].normalBuffer[2] += fnorm[2]*angle;
								}

								{
									a1t[0]=vA[0]-vC[0];a1t[1]=vA[1]-vC[1];a1t[2]=vA[2]-vC[2];
									a2t[0]=vB[0]-vC[0];a2t[1]=vB[1]-vC[1];a2t[2]=vB[2]-vC[2];

									let angle = 0;
									if( (a1t[0]*a1t[0]+a1t[1]*a1t[1]+a1t[2]*a1t[2] ) >0.00000001 && 
										(a2t[0]*a2t[0]+a2t[1]*a2t[1]+a2t[2]*a2t[2] ) >0.00000001 )
										angle = 2*Math.acos( clamp((a1t[0]*a2t[0]+a1t[1]*a2t[1]+a1t[2]*a2t[2])/(Math.sqrt(a1t[0]*a1t[0]+a1t[1]*a1t[1]+a1t[2]*a1t[2])*Math.sqrt(a2t[0]*a2t[0]+a2t[1]*a2t[1]+a2t[2]*a2t[2] ) ), 1.0) );
									pointStateHolder[ci].normalBuffer[0] += fnorm[0]*angle;
									pointStateHolder[ci].normalBuffer[1] += fnorm[1]*angle;
									pointStateHolder[ci].normalBuffer[2] += fnorm[2]*angle;
								}
						}
// --^-^-^-^-^-^-- END GENERATE NORMALS --^-^-^-^-^-^--
					}
				}
			}
		}
	}


	for( var z = 0; z < dim2; z++ ) {
	
		let odd = 0;
		let zOdd = z & 1;

		// for all bounday crossed points, generate the faces from the intersection points.
		for( var y = 0; y < dim1; y++ ) {
			for( var x = 0; x < dim0; x++ ) {

				const normOffset = (x + (y*dim0) + z*dim0*dim1)*5;
				const baseOffset = (x + (y*dim0) + z*dim0*dim1)*6;
				let tetSkip = 0;
				if( x >= (dim0-1)) tetSkip |= 1;
				if( y >= (dim1-1)) tetSkip |= 2;
				if( z >= (dim2-1)) tetSkip |= 4;
				const dataOffset = x + (y*dim0) + z*dim1*dim0;
	        		odd = (( x + y ) &1) ^ zOdd;
				for( tet = 0; tet < 5; tet++ ) {
					if( tetMasks[odd][tet] & tetSkip ) continue;

					let f;
					let invert = 0;
					let useFace = 0;

					// this is 'valid combinations' check.
					if( crosses[ baseOffset+edgeToComp[odd][tet][0] ] ) {
						//console.log( `Output: odd:${odd} tet:${tet} x:${x} y:${y} a:${JSON.stringify(a)}` );
						if( crosses[ baseOffset+edgeToComp[odd][tet][1] ] ) {
							if( crosses[ baseOffset+edgeToComp[odd][tet][2] ] ) {
								// lower left tet. // source point is 0
								useFace = 1;
								invert = ( data[dataOffset+vertToData[odd][tet][0]] >= 0 )?1:0;
							} else {
								if( crosses[ baseOffset+edgeToComp[odd][tet][4] ] && crosses[ baseOffset+edgeToComp[odd][tet][5] ]) {
									// source point is 2? 1?   (0?3?)
									useFace = 2;
									invert = ( data[dataOffset+vertToData[odd][tet][0]] >= 0 )?1:0 ;
								}
							}
						} else {
							if( crosses[ baseOffset+edgeToComp[odd][tet][2]] && crosses[ baseOffset+edgeToComp[odd][tet][3]] && crosses[ baseOffset+edgeToComp[odd][tet][4] ] ) {
								// source point is ? 1? 3?   (0? 2?)
								useFace = 3;
								invert = ( data[dataOffset+vertToData[odd][tet][0]] >= 0 )?1:0  ;
							}else if( crosses[ baseOffset+edgeToComp[odd][tet][3]] && crosses[ baseOffset+edgeToComp[odd][tet][5] ] ) {
								// source point is 1
								useFace = 4;
								invert = ( data[dataOffset+vertToData[odd][tet][1]] >= 0 )?1:0
							}
						}
					} else {
						if( crosses[ baseOffset+edgeToComp[odd][tet][1] ] ) {
							if( crosses[ baseOffset+edgeToComp[odd][tet][2] ] && crosses[ baseOffset+edgeToComp[odd][tet][3] ] && crosses[ baseOffset+edgeToComp[odd][tet][5] ]) {
								// 0?1?   2?3?
								useFace = 5;
								invert = ( data[dataOffset+vertToData[odd][tet][0]] >= 0 )  ?1:0
							} else if( crosses[ baseOffset+edgeToComp[odd][tet][3]] && crosses[ baseOffset+edgeToComp[odd][tet][4] ] ) {
								// source point is 2
								useFace = 6;
								invert = ( data[dataOffset+vertToData[odd][tet][2]] >= 0 ) ?1:0
							}
						} else {
							if( crosses[ baseOffset+edgeToComp[odd][tet][2] ] && crosses[ baseOffset+edgeToComp[odd][tet][4]] && crosses[ baseOffset+edgeToComp[odd][tet][5] ] ) {
								// source point is 3
								useFace = 7;
								invert = ( data[dataOffset+vertToData[odd][tet][3]] >= 0 ) ?1:0
							} else {
							}
						}
					}

					//if( useFace > 5 || useFace < 5 ) continue;
					if( useFace-- ) {
						const fpi = facePointIndexes[odd][tet][invert][useFace];
						if( fpi.length === 1 ){
							const ai = points[baseOffset+fpi[0][0]];
							const bi = points[baseOffset+fpi[0][1]];
							const ci = points[baseOffset+fpi[0][2]];

		
							//console.log( "vertices", tet, useFace, tri, "odd:",odd, "invert:", invert, "pos:", x, y, z, "dels:", pointStateHolder[ai].typeDelta, pointStateHolder[bi].typeDelta, pointStateHolder[ci].typeDelta, "a:", pointStateHolder[ai].invert, pointStateHolder[ai].type1, pointStateHolder[ai].type2, "b:", pointStateHolder[bi].invert, pointStateHolder[bi].type1, pointStateHolder[bi].type2, "c:", pointStateHolder[ci].invert, pointStateHolder[ci].type1, pointStateHolder[ci].type2 );
							const p = [0,0,0], n = [0,0,0];
							TetVert( p, n, pointStateHolder[ai], pointStateHolder[bi], pointStateHolder[ci] );
							normals[normOffset+tet] = {id:0,p:p,n:n, i:invert};
							//console.log( "Setting normal:", x, y, z, dataOffset, normOffset, tet, bits[dataOffset] );
				                }else {
							const ai = points[baseOffset+fpi[0][0]];
							const bi = points[baseOffset+fpi[0][1]];
							const ci = points[baseOffset+fpi[0][2]];

		
							//console.log( "vertices", tet, useFace, tri, "odd:",odd, "invert:", invert, "pos:", x, y, z, "dels:", pointStateHolder[ai].typeDelta, pointStateHolder[bi].typeDelta, pointStateHolder[ci].typeDelta, "a:", pointStateHolder[ai].invert, pointStateHolder[ai].type1, pointStateHolder[ai].type2, "b:", pointStateHolder[bi].invert, pointStateHolder[bi].type1, pointStateHolder[bi].type2, "c:", pointStateHolder[ci].invert, pointStateHolder[ci].type1, pointStateHolder[ci].type2 );
							const p = [0,0,0], n = [0,0,0];

							TetVert( p, n, pointStateHolder[ai], pointStateHolder[bi], pointStateHolder[ci], invert );

							const ai2 = points[baseOffset+fpi[0][0]];
							const bi2 = points[baseOffset+fpi[0][1]];
							const ci2 = points[baseOffset+fpi[0][2]];

		
							//console.log( "vertices", tet, useFace, tri, "odd:",odd, "invert:", invert, "pos:", x, y, z, "dels:", pointStateHolder[ai].typeDelta, pointStateHolder[bi].typeDelta, pointStateHolder[ci].typeDelta, "a:", pointStateHolder[ai].invert, pointStateHolder[ai].type1, pointStateHolder[ai].type2, "b:", pointStateHolder[bi].invert, pointStateHolder[bi].type1, pointStateHolder[bi].type2, "c:", pointStateHolder[ci].invert, pointStateHolder[ci].type1, pointStateHolder[ci].type2 );
							const p2 = [0,0,0], n2 = [0,0,0];
							TetVert( p2, n2, pointStateHolder[ai2], pointStateHolder[bi2], pointStateHolder[ci2] );
	
							p[0] = (p[0]+p2[0])/2;
							p[1] = (p[1]+p2[1])/2;
							p[2] = (p[2]+p2[2])/2;
	                                        	// these are probably going to be pretty close?
							n[0] = (n[0]+n2[0])/2;
							n[1] = (n[1]+n2[1])/2;
							n[2] = (n[2]+n2[2])/2;
							//console.log( "Setting normal:", x, y, z, dataOffset, normOffset, tet, bits[dataOffset] )
							normals[normOffset+tet] = {id:0,p:p,n:n, i:invert};
						}

						bits[dataOffset] = 1;
						content[normOffset + tet] = 1;
					}
				}
			}
		}
	}

	//const dirNums = [ 2,1,7,4,2,4]

	// 0 = Y X Z
	// 1 = X Z Y
	// 2 = Z X Y
	// 3 = x y z
	// 4 = y z x
	// 5 = z y x

	// X 1 3
	// y 0 4
	// z 2 5
	// y > x 0 4 5	
	// x > y 1 2 3 
	// y > z 0 3 4
	// x > z  0 1 3
	// z < x  2 4 5
	// x <  4 5


	function getNormalState( s, baseNormal ) {
		let bnx = baseNormal[0];
		let bny = baseNormal[1];
		let bnz = baseNormal[2];
		if( bnx < 0 ) { s.dir |= 1; bnx = -bnx; }
		if( bny < 0 ) { s.dir |= 2; bny = -bny; }
		if( bnz < 0 ) { s.dir |= 4; bnz = -bnz; }
		if( bnx > bny )  // x > y
			if( bnx > bnz ) // x > z
				if( bny > bnz ) s.largest=0;//0*3+2; // y>z xmax  z min
				else s.largest=1;//0*3+1; // xmax, ymin
			else s.largest=2;//2*3+1; // z longest, x mid  y least
		else if( bny > bnz ) // x < y, y > z   (y > x ) 
			if( bnx>bnz ) s.largest=3;//1*3+2;  // x>z   z min
			else s.largest=4;//1*3+0;  // y longest, z mid, x short
		else s.largest = 5;//2*3 + 0; // z longest, y mid,  x min
	}




	// so at this point I've computed all the cross points
	// all the normals
	// all the centroid/offset points back to the dual phase
	// created the points per-tet-per-cell 
	// now - to actually construct the faces from them.

	function followTheYellowBrickRoad() {
		var here;
		while( here = list.pop() ) {
			for( var nearTet = 0; nearTet < 4; nearTet++ ) ;
		}
	}

	function addFace( ai, bi, ci, n ) {
		if( !n ){
			const tmp1 = [vertices[bi].x-vertices[ai].x, vertices[bi].y-vertices[ai].y, vertices[bi].z-vertices[ai].z];
			const tmp2 = [vertices[ci].x-vertices[ai].x, vertices[ci].y-vertices[ai].y, vertices[ci].z-vertices[ai].z];
			n = cross( [0,0,0], tmp1, tmp2);
			const s = 1/Math.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
			n[0] *= s;
			n[1] *= s;
			n[2] *= s;
		}
		return faces.push( f = new THREE.Face3( ai, bi, ci, new THREE.Vector3(n[0],n[1],n[2]) ) );

	}
	function addPoint( p ) {
		if( !p.id )
			p.id = (vertices.push(new THREE.Vector3( p.p[0],  p.p[1], p.p[2] )),vertices.length-1)+1;

		return p.id - 1;

		return opts.geometryHelper.addPoint( p.p
			 , null, null // texture, and uv[1,0] 
			 , [0xA0,0x00,0xA0,255] // edge color
			 , [0x11, 0x11, 0x11, 255] // face color
			 , [0,0,0] // normal *needs to be updated*;
			 , 100 // pow
			 , false // use texture
			 , false // flat
			 , false // decal texture?
			 , p.p  // modulous of this point.
			 , elements.data[data1], elements.data[data0], t, true
		);
	}
	//https://en.wikipedia.org/wiki/Skew_lines#Nearest_Points

	function getFold( faceNormal, p1, p2, p3, p4 ) {
		p1 = vertices[p1];
		p2 = vertices[p2];
		p3 = vertices[p3];
		p4 = vertices[p4];

		const del1 = [p2.x-p1.x,p2.y-p1.y,p2.z-p1.z];

		const d1 = [p1.x-p3.x,p1.y-p3.y,p1.z-p3.z];
		const d2 = [p2.x-p4.x,p2.y-p4.y,p2.z-p4.z];
		
		const l1 = 1/Math.sqrt(d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2]);
		const l2 = 1/Math.sqrt( d2[0]*d2[0]+d2[1]*d2[1]+d2[2]*d2[2]);
		const n1 = [d1[0]*l1,d1[1]*l1,d1[2]*l1];
		const n2 = [d2[0]*l2,d2[1]*l2,d2[2]*l2];
		const c = [0,0,0];
		cross(c,n1,n2);
		const c2 = [0,0,0];
		cross(c2,n1,c);
	
		const dot1 = del1[0]*c2[0] + del1[1]*c2[1] + del1[2]*c2[2] ;
		const dot2 = n1[0]*c2[0] + n1[1]*c2[1] + n1[2]*c2[2];

		const closeP1 = [p1.x + (( dot1/dot2  )*n1[0]),p1.y + (( dot1/dot2  )*n1[1]),p1.z + (( dot1/dot2  )*n1[2])]

		const dot3 = -del1[0]*c[0] + -del1[1]*c[1] + -del1[2]*c[2] ;
		const dot4 = n2[0]*c[0] + n2[1]*c[1] + n2[2]*c[2];

		const closeP2 = [p2.x + (( dot3/dot4  )*n2[0]),p2.y + (( dot3/dot4  )*n2[1]),p2.z + (( dot3/dot4  )*n2[2])]

		const crossDel = [closeP1[0]-closeP2[0],closeP1[1]-closeP2[1],closeP1[2]-closeP2[2]];

		if( ( crossDel[0] * faceNormal[0] + crossDel[1] * faceNormal[1] + crossDel[2] * faceNormal[2] ) > 0 )
			return 1;
		return 0;

		// compare the length, and prefer to fold on the shorter one.
		if( d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2] < d2[0]*d2[0]+d2[1]*d2[1]+d2[2]*d2[2] ) {
			//console.log( "YES")
			return 1;
		}
		//console.log( "NO", d1, d2)
		return 0;
	}

	if( 0 ) {
		// draw a 'unit' cube for orientation.
	var c1 = addPoint( [0,0,0]);
	var c2 = addPoint( [1,0,0]);
	var c3 = addPoint( [0,2,0]);
	var c4 = addPoint( [1,2,0]);
	var c5 = addPoint( [0,0,4]);
	var c6 = addPoint( [1,0,4]);
	var c7 = addPoint( [0,2,4]);
	var c8 = addPoint( [1,2,4]);
	addFace( c1, c2, c4 );
	addFace( c1, c4, c3 );

	addFace( c2, c6, c8 );
	addFace( c2, c8, c4 );

	addFace( c3, c4, c8 );
	addFace( c3, c8, c7 );
	}

	for( var z = 0; z < dim2; z++ ) {
		//if( z > 3 ) continue;
		let odd = 0;
		let zOdd = z & 1;

		// for all bounday crossed points, generate the faces from the intersection points.
		for( var y = 0; y < dim1; y++ ) {
			//if( y > 3 ) continue;
			for( var x = 0; x < dim0; x++ ) {


				
				const odd = (x+y+z)&1;
				const offset = (x + (y*dim0) + z*dim0*dim1);
				//console.log( "scanning normal:", x, y, z, bits[offset], offset )
				if( !bits[offset] ) {
					// if none of these have any content... 
					//console.log( "Skipping ", x, y, z, offset )
					continue
				}
				let added = 0;
				const baseOffset = offset*5;
				let baseNormal;
				let baseNormalN;

				let n0,n1,n2,n3;
				if( 1  && ( baseNormal = normals[baseOffset+4] ) ) {

					const normDir = { dir:0, largest:0 };
					let inv_ = (normDir.dir&7)?1:0;//normals[baseOffset+4].i;
					let inv = inv_;

					getNormalState( normDir, baseNormal.n );
					// goes through the center.
				if(1)
					if( !content[baseOffset+0] ) {
						// sanity check that 1 2 and 3 are good?
						if(1)
						if( !odd ) {
							if( ( !content[baseOffset+1] )||( !content[baseOffset+2] )||( !content[baseOffset+3] ) )
							{
								console.log( "Bad computation, missing a required intersection if the surface is here..." );
								continue;
							}
							added++;
							switch( normDir.dir ) {
								case 4: case 5: case 6: case 0:
									inv = 0;
									break;
								default:
									inv = 1;
									
							}

							const p0 = addPoint( normals[baseOffset+4]);
							const p1 = addPoint( normals[baseOffset+1]);
							const p2 = addPoint( normals[baseOffset+2]);
							const p3 = addPoint( normals[baseOffset+3]);
							if( inv ){
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
								addFace( p0, p3, p1 );
		
							}else {
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
								addFace( p0, p1, p3 );
							}
						} else {
							if( ( !content[baseOffset+1] )||( !content[baseOffset+2] )||( !content[baseOffset+3] ) )
							{
								console.log( "Bad computation, missing a required intersection if the surface is here..." );
								continue;
							}
							added++;

							switch( normDir.dir ){
								case 0: case 1: case 2: case 4:
									inv = 1;
									break;
								default:
									inv = 0;
							}

							const p0 = addPoint( normals[baseOffset+4]);
							const p1 = addPoint( normals[baseOffset+1]);
							const p2 = addPoint( normals[baseOffset+2]);
							const p3 = addPoint( normals[baseOffset+3]);
							if( inv ){
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
								addFace( p0, p3, p1 );
		
							}else {
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
								addFace( p0, p1, p3 );
							}
						}
					}
					else if( !content[baseOffset+1] ) {
						if(1)
						if( !odd ) {
							if( ( !content[baseOffset+0] )||( !content[baseOffset+2] )||( !content[baseOffset+3] ) )
							{
								console.log( "Bad computation, missing a required intersection if the surface is here..." );
								continue;
							}
							added++;
							switch( normDir.dir ){
								case 1: case 0: case 5: case 3:
									inv = 0;
									break;
								default:
									inv = 1;
							}

							const p0 = addPoint( normals[baseOffset+4]);
							const p1 = addPoint( normals[baseOffset+0]);
							const p2 = addPoint( normals[baseOffset+2]);
							const p3 = addPoint( normals[baseOffset+3]);
							if( inv ) {
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
								addFace( p0, p1, p3 );
							}else{
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
								addFace( p0, p3, p1 );
							}
						} else {
							if( ( !content[baseOffset+0] )||( !content[baseOffset+2] )||( !content[baseOffset+3] ) )
							{
								console.log( "Bad computation, missing a required intersection if the surface is here..." );
								continue;
							}
							added++;
							switch( normDir.dir ){
								case 0: case 2: case 3: case 6:
									inv = 0;
									break;;
								default:
									inv = 1;
									break;
							}

							const p0 = addPoint( normals[baseOffset+4]);
							const p1 = addPoint( normals[baseOffset+0]);
							const p2 = addPoint( normals[baseOffset+2]);
							const p3 = addPoint( normals[baseOffset+3]);
							if( inv ) {
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
								addFace( p0, p1, p3 );
							}else{
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
								addFace( p0, p3, p1 );
							}

						}
						inv = inv_;
					}
					else if( !content[baseOffset+2] ) {
						if(1)
						if( !odd ) {
							if( ( !content[baseOffset+1] )||( !content[baseOffset+0] )||( !content[baseOffset+3] ) )
							{
								console.log( "Bad computation, missing a required intersection if the surface is here..." );
								continue;
							}
							added++;
							inv = 0;
							switch(normDir.dir ){
								case 5:case 4:case 7:case 1:
									break;
								default:
									inv = !inv;
							}							
							const p0 = addPoint( normals[baseOffset+4]);
							const p1 = addPoint( normals[baseOffset+0]);
							const p2 = addPoint( normals[baseOffset+1]);
							const p3 = addPoint( normals[baseOffset+3]);
							if( inv ) {
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
								addFace( p0, p1, p3 );
							}else {
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
								addFace( p0, p3, p1 );
		
							}
						}else {
							if( ( !content[baseOffset+1] )||( !content[baseOffset+0] )||( !content[baseOffset+3] ) )
							{
								console.log( "Bad computation, missing a required intersection if the surface is here..." );
								continue;
							}
							inv = 0;
							switch(normDir.dir ){
								case 6:case 4:case 7:case 2:
									break;
								default:
									inv = !inv;
							}							
							const p0 = addPoint( normals[baseOffset+4]);
							const p1 = addPoint( normals[baseOffset+0]);
							const p2 = addPoint( normals[baseOffset+1]);
							const p3 = addPoint( normals[baseOffset+3]);
							if( inv ) {
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
								addFace( p0, p1, p3 );
							}else {
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
								addFace( p0, p3, p1 );
		
							}
						}
						inv = inv_;
					}
					else if( !content[baseOffset+3] ) {
						if(1)
						if( !odd ) {
							if( ( !content[baseOffset+1] )||( !content[baseOffset+2] )||( !content[baseOffset+0] ) )
							{
								console.log( "Bad computation, missing a required intersection if the surface is here..." );
								continue;
							}
							added++;
							inv = 0;
							switch( normDir.dir) {
								case 0: case 2: case 1: case 4: 
									inv = !inv;
									break;
								default:
									//if( !(normDir.dir & 4) )inv = !inv;
									break;
							}

							const p0 = addPoint( normals[baseOffset+4]);
							const p1 = addPoint( normals[baseOffset+0]);
							const p2 = addPoint( normals[baseOffset+1]);
							const p3 = addPoint( normals[baseOffset+2]);
							if( inv ) {
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
								addFace( p0, p1, p3 );
							}else {
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
								addFace( p0, p3, p1 );
							}


						} else {
							if( ( !content[baseOffset+1] )||( !content[baseOffset+2] )||( !content[baseOffset+0] ) )
							{
								console.log( "Bad computation, missing a required intersection if the surface is here..." );
								continue;
							}
							added++;
							inv = 0;
							switch( normDir.dir) {
								case 3: case 2: case 1: case 7: 
									inv = !inv;
									break;
								default:
									//if( !(normDir.dir & 4) )inv = !inv;
									break;
							}
							const p0 = addPoint( normals[baseOffset+4]);
							const p1 = addPoint( normals[baseOffset+0]);
							const p2 = addPoint( normals[baseOffset+1]);
							const p3 = addPoint( normals[baseOffset+2]);
							if( inv ) {
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
								addFace( p0, p1, p3 );
							}else {
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
								addFace( p0, p3, p1 );
							}

						}
						inv = inv_;
					}
					else {
						// all 5 gets - solid quad emit.
						added++;
						{
							const l = normDir.largest;
							inv = 0;
							const p0 = addPoint( normals[baseOffset+4]);
							const p1 = addPoint( normals[baseOffset+tetCentroidFacet[odd][inv][normDir.dir][l][0]]);
							const p2 = addPoint( normals[baseOffset+tetCentroidFacet[odd][inv][normDir.dir][l][1]]);
							const p3 = addPoint( normals[baseOffset+tetCentroidFacet[odd][inv][normDir.dir][l][2]]);
							const p4 = addPoint( normals[baseOffset+tetCentroidFacet[odd][inv][normDir.dir][l][3]]);
							addFace( p0, p1, p2 );
							addFace( p0, p2, p3 );
							addFace( p0, p3, p4 );
							addFace( p0, p4, p1 );
						}
					}


				if(1)
					if( !odd ) {
						if(1)
						if( (n0 = normals[baseOffset+3]  )) {
							if( (n1 = normals[baseOffset+(1*dim0)*tetCount+1])
							   && (n2 = normals[baseOffset+(1+1*dim0)*tetCount+0] )
							   && (n3 = normals[baseOffset+( 1 )*tetCount+2]) ) {
								added++;
								
								inv = 0;
								const l1 =1/Math.sqrt(  n0.n[0]*n0.n[0] + n0.n[1]*n0.n[1] + n0.n[2]*n0.n[2] ) ;
								const c = [n0.n[0]*l1,n0.n[1]*l1,n0.n[2]*l1];

								const l2 =1/Math.sqrt( n2.n[0]*n2.n[0]+n2.n[1]*n2.n[1]+n2.n[2]*n2.n[2] ) ;
								const c2 = [n2.n[0]*l2,n2.n[1]*l2,n2.n[2]*l2];
								//ncross( c, n2.n, n0.n );
								const d = c[0]*c2[0]+c[1]*c2[1]+c[2]*c2[2];
								if( d > 0.5 ) {
									if( normDir.dir & 4 ) inv = 0; else inv = 1;
									switch( normDir.largest )		{
										//case 0: case 3: case 4: break; // y > z
										//case 0: case 4: case 5: break; // y > x
										case 0: case 1: case 3: {
											if( ( normDir.largest ===3 )||( normDir.largest ===0 ) )
												if( normDir.dir === 0 ) inv = !inv;
											break; // x > z
										}
										default:break;
									}
									const p2 = addPoint( normals[baseOffset+3]);
									const p1 = addPoint( normals[baseOffset + (1*dim0)*tetCount+1]);
									const p0 = addPoint( normals[baseOffset + (1+1*dim0)*tetCount + 0]);
									const p3 = addPoint( normals[baseOffset + ( 1 )*tetCount + 2]);
									
								if( getFold( n0.n, p0, p1, p2, p3 ) )
									if( inv ){
										addFace( p0, p2, p1 );
										addFace( p0, p3, p2 );
									}else {
										addFace( p0, p1, p2 );
										addFace( p0, p2, p3 );
									}
								else
									if( inv ){
										addFace( p1, p3, p2 );
										addFace( p1, p0, p3 );
									}else {
										addFace( p1, p2, p3 );
										addFace( p1, p3, p0 );
									}
								}
							}
						
							if(1)
							if( (n1=normals[baseOffset+0] )
							   && normals[baseOffset + ( 1*dim0*dim1 )*tetCount + 0] 
							   && normals[baseOffset + ( 1*dim0*dim1 )*tetCount + 3] 
							   && normals[baseOffset + ( 1*dim0*dim1 )*tetCount + 4] 
							   ) {  // and 3 because is 4
								added++;
								inv = 0;
								switch( normDir.largest ){
									case 5: case 3: case 4:
										if(( normDir.dir & 2 )) inv = !inv;
										break;
									default:
										if(!( normDir.dir & 1 )) inv = !inv;
										break;
								}

								const p0 = addPoint( normals[baseOffset+0]);
								const p1 = addPoint( normals[baseOffset+3]);
								const p2 = addPoint( normals[baseOffset + ( 1*dim0*dim1)*tetCount + 3]);
								const p3 = addPoint( normals[baseOffset + ( 1*dim0*dim1)*tetCount + 0]);

								if( getFold( n0.n, p0, p1, p2, p3 ) )

									if( inv ){
										addFace( p0, p2, p1 );
										addFace( p0, p3, p2 );
									}else {
										addFace( p0, p1, p2 );
										addFace( p0, p2, p3 );
									}
								else
									if( inv ){
										addFace( p1, p3, p2 );
										addFace( p1, p0, p3 );
									}else {
										addFace( p1, p2, p3 );
										addFace( p1, p3, p0 );
									}
							}

							if(1) // forward in-plane
							if( content[baseOffset + 2] 
									&& normals[baseOffset + (1*dim0)*tetCount + 0]
									&& normals[baseOffset + (1*dim0)*tetCount + 1]
									&& normals[baseOffset + (1*dim0)*tetCount + 4]
									 ) {
								// nothing above, must be... this.
								added++;
								if( normDir.dir & 1 ) inv = 0; else inv = 1;
								switch( normDir.largest ){
									case 0: case 1: case 3:
										break;
									default:
										if( normDir.dir & 1 ){
											if(( normDir.dir & 4 )) inv = !inv;
										}
										else
											if( !(normDir.dir & 4) ) inv = !inv;
										break;

								}
								
								const p0 = addPoint( n0 );
								const p1 = addPoint( normals[baseOffset+2]);
								const p2 = addPoint( normals[baseOffset + (1*dim0)*tetCount + 0]);
								const p3 = addPoint( normals[baseOffset + (1*dim0)*tetCount + 1]);

								if( getFold( n0.n, p0, p1, p2, p3 ) )
									if( inv ){
										addFace( p0, p2, p1 );
										addFace( p0, p3, p2 );
									}else {
										addFace( p0, p1, p2 );
										addFace( p0, p2, p3 );
									}
								else
									if( inv ){
										addFace( p1, p3, p2 );
										addFace( p1, p0, p3 );
									}else {
										addFace( p1, p2, p3 );
										addFace( p1, p3, p0 );
									}

							}

							if(1)
							if( content[baseOffset+1] 
							  && normals[baseOffset + ( 1 )*tetCount + 2] 
							  && normals[baseOffset + ( 1 )*tetCount + 0] 
							  && normals[baseOffset + ( 1 )*tetCount + 4] 
							) {  // and 3 because is 4
								added++;
								inv = 0;
								if( normDir.dir & 2 ) inv = 0; else inv = 1;
								switch( normDir.largest ){
									case 0: case 4: case 3: break;
									default:
										if( normDir.dir & 2 ) {
											if( normDir.dir & 4 )	inv = !inv;
										}else {
											if(!( normDir.dir & 4))	inv = !inv;

										}
								
								}
								const p0 = addPoint( n0 = normals[baseOffset+1]);
								const p1 = addPoint( normals[baseOffset+3]);
								const p2 = addPoint( normals[baseOffset + (1 )*tetCount + 2]);
								const p3 = addPoint( normals[baseOffset + ( 1)*tetCount + 0]);

								if( getFold( n0.n, p0, p1, p2, p3 ) )

									if( inv ){
										addFace( p0, p2, p1 );
										addFace( p0, p3, p2 );
									}else {
										addFace( p0, p1, p2 );
										addFace( p0, p2, p3 );
									}
								else
									if( inv ){
										addFace( p1, p3, p2 );
										addFace( p1, p0, p3 );
									}else {
										addFace( p1, p2, p3 );
										addFace( p1, p3, p0 );
									}

							}

							inv = inv_;
							if(1) // forward vertical
								if( (n0 = normals[baseOffset+3] )
								&& (n1 = normals[baseOffset + (1*dim0)*tetCount+1] )
								&& (n2 = normals[baseOffset + (1*dim0+1*dim0*dim1)*tetCount + 1])
								&& (n3 = normals[baseOffset + (1*dim0*dim1)*tetCount + 3] )) {
									const p0 = addPoint( n0);
									const p1 = addPoint( n1);
									const p2 = addPoint( n2);
									const p3 = addPoint( n3);
									inv = 0;

									switch( normDir.largest ){
										case 0:case 1: case 3:
											if( !(normDir.dir & 1) ) inv = 1;
											break;
										default:
											if( (normDir.dir & 4) ) inv = 1;
											if(! ( (normDir.dir != 2) && (normDir.dir != 5) ) ){
												if( (normDir.dir & 2 ))
													inv = 1;
												else 
													inv = 0;
											};
											break;
									}
									if( getFold( n0.n, p0, p1, p2, p3 ) )
									if( inv ){
										addFace( p0, p2, p1 );
										addFace( p0, p3, p2 );
									}else {
										addFace( p0, p1, p2 );
										addFace( p0, p2, p3 );
									}
								else
									if( inv ){
										addFace( p1, p3, p2 );
										addFace( p1, p0, p3 );
									}else {
										addFace( p1, p2, p3 );
										addFace( p1, p3, p0 );
									}
								}

							if(1) // forawrd-right plane... 
							if( content[baseOffset+ ( 1*dim0*dim1)*tetCount + 3] 
								&& content[baseOffset+ ( 1 + 1*dim0*dim1)*tetCount + 2] 
								&& content[baseOffset + ( 1 )*tetCount + 2] 
								) {
								// right up...
								inv = 0;
								if( (normDir.dir & 2) ) inv = 0; else inv = 1;

								added++;
								
								
								const p0 = addPoint( n0 );
								const p1 = addPoint( normals[baseOffset + (1*dim0*dim1)*tetCount + 3]);
								const p2 = addPoint( normals[baseOffset + (1+1*dim0*dim1)*tetCount + 2]);
								const p3 = addPoint( normals[baseOffset + ( 1 )*tetCount + 2]);
								if( getFold( n0.n, p0, p1, p2, p3 ) )
									if( inv ){
										addFace( p0, p2, p1 );
										addFace( p0, p3, p2 );
									}else {
										addFace( p0, p1, p2 );
										addFace( p0, p2, p3 );
									}
								else
									if( inv ){
										addFace( p1, p3, p2 );
										addFace( p1, p0, p3 );
									}else {
										addFace( p1, p2, p3 );
										addFace( p1, p3, p0 );
									}
							}


						}
						// to the fore
						// TODO : THIS IS STILL BROKEN...

///----------------------------------------------------------
//    ODD  ....   (above is even (!odd) )
///----------------------------------------------------------

					}else {
						// all other cases that don't include 4 will be rendered in another cell...
						// this one can only add new faces here.
						// to the right

						// directly up.
						if(1)
						if( content[baseOffset+1] 
							&& content[baseOffset+2] 
						   && normals[baseOffset + ( 1*dim0*dim1 )*tetCount + 2] 
						   && normals[baseOffset + ( 1*dim0*dim1 )*tetCount + 4] 
						   && normals[baseOffset + ( 1*dim0*dim1 )*tetCount + 1] 
						  ) {  // and 3 because is 4
							added++;
							inv = 0;
							//if( normDir.dir & 1 ) inv = 0; else inv = 1;
							switch( normDir.largest ){
								//case 1:case 3: case 0:
								//case 4:case 3: case 0:
								case 4: case 3: case 5:
									if(!( normDir.dir & 2) ) inv = !inv;
									break;
									default:
										if(!( normDir.dir & 1) ) inv = !inv;
										break;
							}

							const p0 = addPoint( n0=normals[baseOffset+1]);
							const p1 = addPoint( normals[baseOffset+2]);
							const p2 = addPoint( normals[baseOffset + ( 1*dim0*dim1)*tetCount + 2]);
							const p3 = addPoint( normals[baseOffset + ( 1*dim0*dim1)*tetCount + 1]);
								if( getFold( n0.n, p0, p1, p2, p3 ) )
									if( inv ){
										addFace( p0, p2, p1 );
										addFace( p0, p3, p2 );
									}else {
										addFace( p0, p1, p2 );
										addFace( p0, p2, p3 );
									}
								else
									if( inv ){
										addFace( p1, p3, p2 );
										addFace( p1, p0, p3 );
									}else {
										addFace( p1, p2, p3 );
										addFace( p1, p3, p0 );
									}
							}


	
						inv = inv_;
						if(1)
						if( n0 = normals[baseOffset+3] ) { 
							// and content[baseOffset+(1*dim0)*tetCount+0]  content[baseOffset+(1*dim0)*tetCount+3] content[baseOffset+(1*dim0)*tetCount+4]
							//console.log( "Content Here:", odd, inv, content[baseOffset+0], content[baseOffset+1], content[baseOffset+2], content[baseOffset+3],content[baseOffset+4] );
							inv = inv_;
							// to the right - in-plane
							if(1)
							if( content[baseOffset+(1*dim0)*tetCount+1] 
							   && content[baseOffset+(1+1*dim0)*tetCount+0] 
							   && content[baseOffset+( 1 )*tetCount+2] ) {
								added++;
								inv = 0;
								//console.log( "SOmething:",content[baseOffset+ ( 1*dim0)*tetCount + 0] ,content[baseOffset+ ( 1*dim0*dim1)*tetCount +1] ,content[baseOffset+ ( 1*dim0*dim1)*tetCount + 2] ,content[baseOffset+ ( 1*dim0*dim1)*tetCount + 3] ,content[baseOffset+ ( 1*dim0*dim1)*tetCount + 4]  )
								if( normDir.dir & 4 ) inv = !inv;
	
								const p0 = addPoint( n0 = normals[baseOffset+3]);
								const p1 = addPoint( normals[baseOffset + ( 1 )*tetCount + 2]);
								const p2 = addPoint( normals[baseOffset + (1+1*dim0)*tetCount + 0]);
								const p3 = addPoint( normals[baseOffset + (1*dim0)*tetCount+1]);
								inv = !inv;
								if( getFold( n0.n, p0, p1, p2, p3 ) )
									if( inv ){
										addFace( p0, p2, p1 );
										addFace( p0, p3, p2 );
									}else {
										addFace( p0, p1, p2 );
										addFace( p0, p2, p3 );
									}
								else
									if( inv ){
										addFace( p1, p3, p2 );
										addFace( p1, p0, p3 );
									}else {
										addFace( p1, p2, p3 );
										addFace( p1, p3, p0 );
									}
							}


							if(1)
							if( content[baseOffset+1]
										&& normals[baseOffset + ( 1 )*tetCount + 2] 
										&& normals[baseOffset + ( 1 )*tetCount + 4] 
										&& normals[baseOffset + ( 1 )*tetCount + 0] 
										) {  // and 3 because is 4
								added++;

								if( normDir.dir & 2 ) inv = 1; else inv = 0;
								switch(normDir.largest  ){
									case 0: case 3: case 4: {
										break;
									};
									default:
										if( (normDir.dir & 2 ) ) {
											if(!( normDir.dir & 4 )) inv = !inv;
										}
										else {
											if(( normDir.dir & 4 )) inv = !inv;
										}
								}
	
								const p0 = addPoint( n0 = normals[baseOffset+1]);
								const p1 = addPoint( normals[baseOffset+3]);
								const p2 = addPoint( normals[baseOffset + ( 1 )*tetCount + 2]);
								const p3 = addPoint( normals[baseOffset + ( 1 )*tetCount + 0]);
								if( getFold( n0.n, p0, p1, p2, p3 ) )
									if( inv ){
										addFace( p0, p2, p1 );
										addFace( p0, p3, p2 );
									}else {
										addFace( p0, p1, p2 );
										addFace( p0, p2, p3 );
									}
								else
									if( inv ){
										addFace( p1, p3, p2 );
										addFace( p1, p0, p3 );
									}else {
										addFace( p1, p2, p3 );
										addFace( p1, p3, p0 );
									}
							}
							// to the fore
								

						}

						// 
						// odd 2 - above/forward line across
						inv = inv_;
						if(1)
						if( n0 = normals[baseOffset + 2] ) {

							if(1)	
							if( (n1 = normals[baseOffset + (1*dim0)*tetCount + 0] )
									&& (n2 =normals[baseOffset + (1*dim0)*tetCount + 1] )
									&& (n3 =normals[baseOffset + 3] )
									&& (normals[baseOffset + (1*dim0)*tetCount + 4] ) // make sure the surface goes here...
									 ) {

								// nothing above, must be... this.
								added++;
								inv = 0;
								//if( !(normDir.dir &4)  ) inv = 1; else inv = 0;
								switch(normDir.dir  ){
									case 7: case 5:
										inv = 1;
										break;;
									case 3: case 1:
										if( normDir.largest === 3 || normDir.largest === 1|| normDir.largest === 0 )
											inv = 1;
										break;
									case 2: case 0:
										break;
									case 6: case 4:
										if( normDir.largest === 2 || normDir.largest === 5 || normDir.largest === 4 ) 
											inv = 1;
										break;
								}

								const p0 = addPoint( n0);
								const p1 = addPoint( n1);
								const p2 = addPoint( n2);
								const p3 = addPoint( n3);
								//console.log( "Emit:",normDir, inv )
								if( getFold( n0.n, p0, p1, p2, p3 ) )
									if( inv ){
										addFace( p0, p2, p1 );
										addFace( p0, p3, p2 );
									}else {
										addFace( p0, p1, p2 );
										addFace( p0, p2, p3 );
									}
								else
									if( inv ){
										addFace( p1, p3, p2 );
										addFace( p1, p0, p3 );
									}else {
										addFace( p1, p2, p3 );
										addFace( p1, p3, p0 );
									}
							}

							if(1)
							if( (n1 = normals[baseOffset+(1*dim0*dim1)*tetCount+2])
								&& (n2 = normals[baseOffset + (1*dim0+1*dim0*dim1)*tetCount + 0])
								//&& ( normals[baseOffset + (1*dim0+1*dim0*dim1)*tetCount + 4])
								&& (n3 = normals[baseOffset + (1*dim0)*tetCount + 0])
								 ) {

								// nothing above, must be... this.
								added++;
								inv  = 0 ;
								switch( normDir.largest ) {
									case 1: case 3: case 2:case 0:
										switch( normDir.dir ) {
											case 5:// case 6: case 7: case 4: 
											case 7: case 1: case 3:
												inv = 1;
												break;
											default:
												//console.log( "This is sometimes inv", x, y, z, odd, inv, normDir );
												//inv = 1;
												break;
										}
										break;
									default:
										if( normDir.dir & 1) {
											inv = 1;
										}
										break;
								}
	
								const p0 = addPoint( n0);
								const p1 = addPoint( n1);
								const p2 = addPoint( n2);
								const p3 = addPoint( n3);
								if( getFold( n0.n, p0, p1, p2, p3 ) )
									if( inv ){
										addFace( p0, p2, p1 );
										addFace( p0, p3, p2 );
									}else {
										addFace( p0, p1, p2 );
										addFace( p0, p2, p3 );
									}
								else
									if( inv ){
										addFace( p1, p3, p2 );
										addFace( p1, p0, p3 );
									}else {
										addFace( p1, p2, p3 );
										addFace( p1, p3, p0 );
									}
							}

							// forward diagonal of 2... 
						}

						inv = inv_;
						if(1) // good!
						// odd, 'up'(y), x, 
						if( n0 = normals[baseOffset+1] ) {
							if( (n1 = normals[baseOffset+ ( 1*dim0*dim1)*tetCount + 1] ) 
								&& (n2 = normals[baseOffset+ ( 1 + 1*dim0*dim1)*tetCount + 0] )
								//&& ( normals[baseOffset+ ( 1 + 1*dim0*dim1)*tetCount + 4] )
								&& normals[baseOffset + ( 1 )*tetCount + 0]
								) {
									const l1 =1/Math.sqrt(  n0.n[0]*n0.n[0] + n0.n[1]*n0.n[1] + n0.n[2]*n0.n[2] ) ;
									const c = [n0.n[0]*l1,n0.n[1]*l1,n0.n[2]*l1];

									const l2 =1/Math.sqrt( n2.n[0]*n2.n[0]+n2.n[1]*n2.n[1]+n2.n[2]*n2.n[2] ) ;
									const c2 = [n2.n[0]*l2,n2.n[1]*l2,n2.n[2]*l2];
									//ncross( c, n2.n, n0.n );
									const d = c[0]*c2[0]+c[1]*c2[1]+c[2]*c2[2];
									//if( c[0]*c[0]+c[1]*c[1]+c[2]*c[2] < 0.5 )
									if( d > 0.8 )
									{
										// right up...
										// make sure this is a complete thing...
										added++;
										if( normDir.dir & 2 ) inv = 0; else inv = 1;

										const p0 = addPoint( n0 = normals[baseOffset+1]);
										const p1 = addPoint( normals[baseOffset + (1*dim0*dim1)*tetCount + 1]);
										const p2 = addPoint( n2);
										const p3 = addPoint( normals[baseOffset + ( 1 )*tetCount + 0]);
										if( getFold( n0.n, p0, p1, p2, p3 ) )
											if( inv ){
												addFace( p0, p2, p1 );
												addFace( p0, p3, p2 );
											}else {
												addFace( p0, p1, p2 );
												addFace( p0, p2, p3 );
											}
										else
											if( inv ){
												addFace( p1, p3, p2 );
												addFace( p1, p0, p3 );
											}else {
												addFace( p1, p2, p3 );
												addFace( p1, p3, p0 );
											}
									}
							}

							
						}
					}
				} else {
					let n0,n1,n2,n3;
				if(1)
					if( !odd ) 
					{
						if(1)
						// not through the center... 
						if( (n2 = normals[baseOffset+3] )
						  && (n1 = normals[baseOffset + (1*dim0)*tetCount + 1] )
						  && (n0 = normals[baseOffset + (1 + 1*dim0)*tetCount + 0])
						  //&& (normals[baseOffset + (1 + 1*dim0)*tetCount + 4])
						  && (n3 = normals[baseOffset + ( 1 )*tetCount + 2] )) {
							  
							const normDir = { dir:0, largest:0 };
							getNormalState( normDir, n0.n );
							let inv = 0;//!(( normDir.dir ^ 0x7 ) & 7);

							switch( normDir.dir ){
								case 0: case 1: case 2: case 3:
										inv = 1;
								break;
								default:
									break;;
							}

						const p0 = addPoint( n0);
							const p1 = addPoint( n1);
							const p2 = addPoint( n2);
							const p3 = addPoint( n3);
							if( getFold( n0.n, p0, p1, p2, p3 ) )
							if( inv ){
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
							}else {
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
							}
						else
							if( inv ){
								addFace( p1, p3, p2 );
								addFace( p1, p0, p3 );
							}else {
								addFace( p1, p2, p3 );
								addFace( p1, p3, p0 );
							}
						}
				//console.log( "status:", x, y, z, baseOffset, n0, n1, n2, n3 );
						if(1)
						if( (n0 = normals[baseOffset+3] )
						  && (n3 = normals[baseOffset + (1*dim0*dim1)*tetCount + 3] )
						  && (n2 = normals[baseOffset + (1*dim0 + 1*dim0*dim1)*tetCount +1])
						  && (n2 = normals[baseOffset + (1*dim0 + 1*dim0*dim1)*tetCount +4])
						  && (n1 = normals[baseOffset + (1*dim0 )*tetCount + 1] )) {
							  // this feels large.

							const normDir = { dir:0, largest:0 };
							getNormalState( normDir, n0.n );
							let inv = 0;//normDir.dir&7;

							switch( normDir.dir ){
								case 0: case 2: case 4: case 6:
										inv = 1;
								break;
								default:
									break;;
							}

							  const p0 = addPoint( n0);
							const p1 = addPoint( n1);
							const p2 = addPoint( n2);
							const p3 = addPoint( n3);
							if( getFold( n0.n, p0, p1, p2, p3 ) )
							if( inv ){
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
							}else {
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
							}
						else
							if( inv ){
								addFace( p1, p3, p2 );
								addFace( p1, p0, p3 );
							}else {
								addFace( p1, p2, p3 );
								addFace( p1, p3, p0 );
							}
						}
				//console.log( "status:", x, y, z, baseOffset, n0, n1, n2, n3 );
						if(1)
						if( (n0 = normals[baseOffset+3] )
						  && (n1 = normals[baseOffset + (1*dim0*dim1)*tetCount + 3] )
						  && (n2 = normals[baseOffset + (1 + 1*dim1*dim0)*tetCount + 2])
						  ///&& ( normals[baseOffset + (1 + 1*dim1*dim0)*tetCount + 4])
						  && (n3 = normals[baseOffset + ( 1 )*tetCount + 2])) {
							  
							  const normDir = { dir:0, largest:0 };
							  getNormalState( normDir, n0.n );
							  let inv = 0;//(normDir.dir&7);
							  switch( normDir.dir ){
								case 0: case 1: case 4: case 5:
										inv = 1;
								break;
								default:
									break;;
							}

							const p0 = addPoint( n0);
							const p1 = addPoint( n1);
							const p2 = addPoint( n2);
							const p3 = addPoint( n3);
							if( getFold( n0.n, p0, p1, p2, p3 ) )
							if( inv ){
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
							}else {
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
							}
						else
							if( inv ){
								addFace( p1, p3, p2 );
								addFace( p1, p0, p3 );
							}else {
								addFace( p1, p2, p3 );
								addFace( p1, p3, p0 );
							}
						}
//				console.log( "status:", x, y, z,baseOffset,  n0, n1, n2, n3 );
					} else {
						// THERE ARE 2 inverted triangels in this case.  TODO
						if(1)
						if( (n0 = normals[baseOffset+1] )
						  && (n1 = normals[baseOffset + (1*dim0*dim1)*tetCount + 1] )
						  && (n2 = normals[baseOffset + (1 + 1*dim0*dim1)*tetCount + 0])
						  && (n3 = normals[baseOffset + ( 1 )*tetCount + 0])) {
							  
							  const normDir = { dir:0, largest:0 };
							  getNormalState( normDir, n0.n );
							  let inv = 1;
							  switch( normDir.dir ){
								case 0: case 1: case 2: case 3:
										inv = 0;
								break;
								default:
									break;;
							}
							  const p0 = addPoint( n0);
							const p1 = addPoint( n1);
							const p2 = addPoint( n2);
							const p3 = addPoint( n3);
							if( getFold( n0.n, p0, p1, p2, p3 ) )
							if( inv ){
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
							}else {
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
							}
						else
							if( inv ){
								addFace( p1, p3, p2 );
								addFace( p1, p0, p3 );
							}else {
								addFace( p1, p2, p3 );
								addFace( p1, p3, p0 );
							}
						}
						// there are two more triangles inverted on this TODO
						if(1)
						if( (n0 = normals[baseOffset+2] )
						  && (n1 = normals[baseOffset + (1*dim0*dim1)*tetCount + 2] )
						  && (n2 = normals[baseOffset + (1*dim0*dim1 + 1*dim0)*tetCount + 0])
						  //&& (normals[baseOffset + (1*dim0*dim1 + 1*dim0)*tetCount + 4])
						  && (n3 = normals[baseOffset + ( 1*dim0 )*tetCount + 0])) {

							const normDir = { dir:0, largest:0 };
							  getNormalState( normDir, n0.n );
							  let inv = 1;//normDir.dir&7;
							  switch( normDir.dir ){
								case 0: case 4: case 2: case 6:
										inv = 0;
								break;
								default:
									break;;
							}
							  const p0 = addPoint( n0);
							const p1 = addPoint( n1);
							const p2 = addPoint( n2);
							const p3 = addPoint( n3);
							if( getFold( n0.n, p0, p1, p2, p3 ) )
							if( inv ){
								addFace( p0, p2, p1 );
								addFace( p0, p3, p2 );
							}else {
								addFace( p0, p1, p2 );
								addFace( p0, p2, p3 );
							}
						else
							if( inv ){
								addFace( p1, p3, p2 );
								addFace( p1, p0, p3 );
							}else {
								addFace( p1, p2, p3 );
								addFace( p1, p3, p0 );
							}
						}
						
						if(1)
						if( (n0 = normals[baseOffset+3] )
						  && (n1 = normals[baseOffset + (1*dim0)*tetCount + 1] )
						  && (n2 = normals[baseOffset + (1*dim0 + 1)*tetCount + 0])
						  //&& (normals[baseOffset + (1*dim0 + 1)*tetCount + 4])
						  && (n3 = normals[baseOffset + ( 1 )*tetCount + 2])) {

							const normDir = { dir:0, largest:0 };
							  getNormalState( normDir, n0.n );
							  let inv = 1;//normDir.dir&7;
							  switch( normDir.dir ){
								case 0: case 1: case 2: case 3:
										inv = 0;
								break;
								default:
									break;;
							}
							const p0 = addPoint( n0);
							const p1 = addPoint( n1);
							const p2 = addPoint( n2);
							const p3 = addPoint( n3);
							if( getFold( n0.n, p0, p1, p2, p3 ) )
								if( inv ){
									addFace( p0, p2, p1 );
									addFace( p0, p3, p2 );
								}else {
									addFace( p0, p1, p2 );
									addFace( p0, p2, p3 );
								}
							else
								if( inv ){
									addFace( p1, p3, p2 );
									addFace( p1, p0, p3 );
								}else {
									addFace( p1, p2, p3 );
									addFace( p1, p3, p0 );
								}
						}
						
						
				//console.log( "status:", x, y, z, n0, n1, n2, n3 );
					}
				}
				//console.log( "ADDED:", added );
			}
		}
	}




	// update geometry (could wait for index.html to do this?
	//if( showGrid )
	//	opts.geometryHelper.markDirty();

	opts.points  = points;   
	opts.normals = normals;
	opts.bits    = bits; 
	// internal utility function to limit angle
	function clamp(a,b) {
		if( a < b ) return a; return b;
	}
	//return vertices;
}

}
})()

if("undefined" != typeof exports) {
	exports.mesher = MarchingTetrahedra4;
}

