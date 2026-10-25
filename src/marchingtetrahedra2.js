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
 *  
 * (Several bug fixes were made to deal with oriented faces)
 *
 * Javascript port by d3x0r
 */
var MarchingTetrahedra2 = (function() {

const _debug = false;

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

const tetTable = [ [
	[ 1,0,5,3 ],
	[ 2,3,6,0 ],
	[ 4,5,0,6 ],
	[ 7,3,5,6 ],
	[ 0,6,5,3 ],
]

,[
	[ 0,2,4,1 ],
	[ 3,1,7,2 ],
	[ 5,7,1,4 ],
	[ 6,7,4,2 ],
	[ 4,7,1,2 ],
] ]

const v = [0,0,0,0];
const g = [null,null,null,null];
const cc = [0,0,0,0]; // cube corner (0-7) of each tet corner, for element data
let cellBase = 0;      // sample index of the cell's (x,y,z-1) corner
let cornerOfs = null;  // sample offset of each cube corner from cellBase

const values = [[-1,-1,-1,-1],[-1,-1,-1,-1]];

const cellOrigin = [0,0,0];

// Per-vertex element data, kept between calls (grown by doubling) so meshing
// doesn't make garbage for it; results are views, valid until the next call.
let mat1 = new Uint16Array(1024), mat2 = new Uint16Array(1024), matDel = new Float32Array(1024);
function growMaterials( n ) {
	const a = new Uint16Array(n), b = new Uint16Array(n), c = new Float32Array(n);
	a.set( mat1 ); b.set( mat2 ); c.set( matDel );
	mat1 = a; mat2 = b; matDel = c;
}


return function(data, dims, opts) {
	opts = opts || { maximize:false,minimize:false,inflation:0};
	const dim0 = dims[0];	
	const dim1 = dims[1];	
	const dim01 = dim0*dim1;
	const dim2 = dims[2];	
   const vertices = []
    , faces = [];
   // per vertex, like MT3/MT4 hand the grid buffer: the elements at the two
   // ends of the vertex's edge, and how far along it the crossing is.
   const elem = opts.elements && opts.elements.data;
   let t1 = 0, t2 = 0, td = 0;


/*
// inverted not applies to that point in or out of shape vs the others.
   0 _ _ 1  (inverted)
   |\   /|
    \\2//     (above page)
    | | |
     \|/
      3  (inverted)
*/


function e2(p) {
	faces.push(p);
}

function emit( p ) {
	const n = vertices.length;
	if( n >= mat1.length ) growMaterials( mat1.length*2 );
	mat1[n] = t1; mat2[n] = t2; matDel[n] = td;
	vertices.push( p );
	return n;
}

// the crossing on the tet edge corner a -> corner b, 'del' of the way along.
// Elements follow MT4: its point is p1+(p0-p1)*t with types (p0, p1, t), so
// here (b, a, del).  emit() picks them up.
function lerp( a, b, del ) {
	const p1 = g[a], p2 = g[b];
	if( elem ) { t1 = elem[cellBase+cornerOfs[cc[b]]]; t2 = elem[cellBase+cornerOfs[cc[a]]]; td = del; }
	return [ cellOrigin[0] + p1[0] + (p2[0]-p1[0])*del
               , cellOrigin[1] + p1[1] + (p2[1]-p1[1])*del
               , cellOrigin[2] + p1[2] + (p2[2]-p1[2])*del ];
}

function tetCompute( values, geometry, invert ) {

	_debug && console.log( "tet: v:", values, "g:", geometry );
	if( ( values[0] >= 0 ) ) {
		invert = !invert;
		values[0] = 0-values[0];
		values[1] = 0-values[1];
		values[2] = 0-values[2];
		values[3] = 0-values[3];
	}

	{
		// 0 is outside
		if( values[1] >= 0 ) {
        		// 1 is inside  0-1 crosses
			const cross1 = ( !( values[1]-values[0] ) )?0: -values[0] / ( values[1]-values[0] );                        
	                if( values[2] >= 0 ) {
                        	// 0-2 is also a cross
        			const cross2 = ( !( values[2]-values[0] ) )?0: -values[0] / ( values[2] - values[0] );
                                if( values[3] >= 0 ) {
		                	// 0-3 is also a cross
        				const cross3 = ( !( values[3]-values[0] ) )?0: -values[0] / ( values[3] - values[0] );
	                                // emit tri.  
                                        if( invert ) {
                                                e2([emit( lerp( 0, 1, cross1 ) ),
		                                    emit( lerp( 0, 2, cross2 ) ),
		                                    emit( lerp( 0, 3, cross3 ) )]);
	                                } else {
                                                e2([emit( lerp( 0, 2, cross2 ) ),
	                                            emit( lerp( 0, 1, cross1 ) ),
	                                            emit( lerp( 0, 3, cross3 ) )] );
	                                }
                        	} else {
        				const cross3 = ( !( values[1]-values[3] ) )?0: -values[3] / ( values[1] - values[3] );
        				const cross4 = ( !( values[2]-values[3] ) )?0: -values[3] / ( values[2] - values[3] );
                                        // emit quad
        				const a= emit(lerp( 0, 1, cross1 ));
					const b= emit(lerp( 0, 2, cross2 ));
					const c= emit(lerp( 3, 1, cross3 ));
					const d= emit(lerp( 3, 2, cross4 ));
	                                // emit a,b,c  b,c,d
                                        if( invert ) {
        					//e2( [a,b,d,c] );
	                                        e2( [c,a,d] );
		                                e2( [d,a,b] );
                                        }else{
        					//e2( [a,c,d,b] );
 		                                e2( [a,c,d] );
                                                e2( [a,d,b] );
					}
	                        }
                        } else {
                        	if( values[3] >= 0 ) {
        				const cross2 = ( !( values[1]-values[2] ) )?0: -values[2] / ( values[1] - values[2] );
					const cross3 = ( !( values[3]-values[0] ) )?0: -values[0] / ( values[3] - values[0] );
					const cross4 = ( !( values[3]-values[2] ) )?0: -values[2] / ( values[3] - values[2] );
	                        	// emit quad
                                        // emit a,b,c  b,c,d
        				const a= emit(lerp( 0, 1, cross1 ));
					const b= emit(lerp( 0, 3, cross3 ));
					const c= emit(lerp( 2, 1, cross2 ));
					const d= emit(lerp( 2, 3, cross4 ));

					// VERIFIED
	                                if( invert ) {
        					//e2( [b,a,c,d] );
	                                        e2( [b,a,c] );
		                                e2( [d,b,c] );
        				}else {
						//e2( [b,d,c,a] );
	        	                        e2( [a,b,c] );
                        	                e2( [b,d,c] );
        				}
	                        } else {
                                	// 0 out  1 in  2 out  3 out
        				if( !( values[1]-values[2] ) )
						cross2 = 0;
					else
	                                        cross2 = -values[2] / ( values[1] - values[2] );
					if( !( values[1]-values[3] ) )
						cross3 = 0;
					else
	                                        cross3 = -values[3] / ( values[1] - values[3] );
	                                // emit tri 2,3,0
	                                if( !invert ) {
                                                 e2([emit(lerp( 2, 1, cross2 )),
		                                         emit(lerp( 0, 1, cross1 )),
		                                         emit(lerp( 3, 1, cross3 ))]);
	                                } else {
                                                 e2([emit(lerp( 0, 1, cross1 )),
	        	                                 emit(lerp( 2, 1, cross2 )),
	                	                         emit(lerp( 3, 1, cross3 ))]);
	                                }
                                }
                        }
                } else {
                	// 0,1 outside
                        if( values[2] >= 0 ) {
                        	// 0-2 is also a cross
        			const cross1 = ( !( values[2]-values[0] ) )?0: -values[0] / ( values[2] - values[0] );
				const cross2 = ( !( values[2]-values[1] ) )?0: -values[1] / ( values[2] - values[1] );
                                if( values[3] >= 0 ) {
		                	// 0-3 is also a cross
        				const cross3 = ( !( values[3]-values[0] ) )?0: -values[0] / ( values[3] - values[0] );
        				const cross4 = ( !( values[3]-values[1] ) )?0: -values[1] / ( values[3] - values[1] );
                                        // emit quad.  
        				const a= emit(lerp( 0, 2, cross1 ));
					const b= emit(lerp( 1, 2, cross2 ));
					const c= emit(lerp( 0, 3, cross3 ));
					const d= emit(lerp( 1, 3, cross4 ));
	                                if( !invert ) {
        					//e2( [a,b,d,c] );
	                                        e2( [d,a,b] );
		                                e2( [c,a,d] );
        				}else {
						//e2( [a,c,d,b] );
	        	                        e2( [a,d,b] );
                        	                e2( [a,c,d] );
        				}
	                	} else {
        				// 0 out 1 out   2 in  3 out
					const cross3 = ( !( values[2]-values[3] ) )?0: -values[3] / ( values[2] - values[3] );
                                        // emit tri 0,1,3
        
                                        if( invert ) {
                                                e2( [ emit( lerp( 1, 2, cross2 ) ),
		                                        emit( lerp( 0, 2, cross1 ) ),
		                                        emit( lerp( 3, 2, cross3 ) ) ] );
	                                } else {
                                                e2( [ emit(  lerp( 0, 2, cross1 ) ),
	                                         emit(lerp( 1, 2, cross2 ) ),
	                                         emit(lerp( 3, 2, cross3 ) )]);
	                                }
                                }
                        } else {
                                // 0,1,2 outside
                        	if( values[3] >= 0 ) {
        				// 3 inside...
					const cross1 = ( !( values[3]-values[0] ) )?0: -values[0] / ( values[3] - values[0] );
        				const cross2 = ( !( values[3]-values[1] ) )?0: -values[1] / ( values[3] - values[1] );
					const cross3 = ( !( values[3]-values[2] ) )?0: -values[2] / ( values[3] - values[2] );
	                        	// emit tri
                                        if( invert ) {
                                                e2( [ emit( lerp( 0, 3, cross1 ) ),
	                                              emit(  lerp( 1, 3, cross2 ) ),
	                                              emit(  lerp( 2, 3, cross3 ) )]);
	                                } else {
                                                e2( [ emit( lerp( 1, 3, cross2 ) ),
	                                              emit(  lerp( 0, 3, cross1 ) ),
	                                              emit(  lerp( 2, 3, cross3 ) )]);
	                                }
                                }
                        }
                }
        } 
}

	// a sample of exactly 0 is outside, like the other meshers and the element
	// data (which is void there); -0 would test as inside below.
	// sample offset of each cube corner from the cell's (x,y,z-1) sample
	cornerOfs = geom.map( c=>c[0] + c[1]*dim0 + c[2]*dim01 );

	// cells are [x,x+1] x [y,y+1] x [z-1,z]: values[0] is the previous z layer.
	// (x and y stop one short -- x+1/y+1 read the next row -- and z=0 only
	// loads a layer; it used to mesh against the previous column's last one.)
	for( let x = 0; x < dim0-1; x++ ) {
		cellOrigin[0] = x;
		for( let y = 0; y < dim1-1; y++ ) {
			cellOrigin[1] = y;
			for( let z = 0; z < dim2; z++ ) {
				cellOrigin[2] = z-1;
				let d;
				const tmp = values[0]; // swap rows
				values[0] = values[1];
				values[1] = tmp;
				
				d = data[ z * dim01 + (y+0) * dim0 + (x+0) ]; values[1][0] = d ? -d : -1e-30;
				d = data[ z * dim01 + (y+0) * dim0 + (x+1) ]; values[1][1] = d ? -d : -1e-30;

				d = data[ z * dim01 + (y+1) * dim0 + (x+0) ]; values[1][2] = d ? -d : -1e-30;
				d = data[ z * dim01 + (y+1) * dim0 + (x+1) ]; values[1][3] = d ? -d : -1e-30;
				if( !z ) continue;

				const alt = (x+y+z)&1;
				cellBase = (z-1) * dim01 + y * dim0 + x;
				for( let tet of tetTable[alt] ) {
					for( let i = 0; i < 4; i++ ) {
						const c = tet[i];
						v[i] = values[c>>2][c&3];
						g[i] = geom[c];
						cc[i] = c;
					}
					tetCompute( v, g, false );
				}
			}
		}
	}



  
  const nv = vertices.length;
  return { vertices: vertices, faces: faces
         , types1: mat1.subarray(0,nv), types2: mat2.subarray(0,nv), deltas: matDel.subarray(0,nv) };
}
})();


if("undefined" != typeof exports) {
  exports.mesher = MarchingTetrahedra2;
}
