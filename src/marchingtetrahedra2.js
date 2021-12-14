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

return function(data, dims, opts) {
	opts = opts || { maximize:false,minimize:false,inflation:0};
   var vertices = []
    , faces = [];


/*
// inverted not applies to that point in or out of shape vs the others.
   0 _ _ 1  (inverted)
   |\   /|
    \\2//     (above page)
    | | |
     \|/
      3  (inverted)
*/

const cellOrigin = [0,0,0];

function e2(p) {
	faces.push(p);
}

function emit( p ) {
	vertices.push( p );
	return vertices.length-1;
}

function lerp( p1, p2, del ) {
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
                                                e2([emit( lerp( geometry[0], geometry[1], cross1 ) ),
		                                    emit( lerp( geometry[0], geometry[2], cross2 ) ),
		                                    emit( lerp( geometry[0], geometry[3], cross3 ) )]);
	                                } else {
                                                e2([emit( lerp( geometry[0], geometry[2], cross2 ) ),
	                                            emit( lerp( geometry[0], geometry[1], cross1 ) ),
	                                            emit( lerp( geometry[0], geometry[3], cross3 ) )] );
	                                }
                        	} else {
        				const cross3 = ( !( values[1]-values[3] ) )?0: -values[3] / ( values[1] - values[3] );
        				const cross4 = ( !( values[2]-values[3] ) )?0: -values[3] / ( values[2] - values[3] );
                                        // emit quad
        				const a= emit(lerp( geometry[0], geometry[1], cross1 ));
					const b= emit(lerp( geometry[0], geometry[2], cross2 ));
					const c= emit(lerp( geometry[3], geometry[1], cross3 ));
					const d= emit(lerp( geometry[3], geometry[2], cross4 ));
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
        				const a= emit(lerp( geometry[0], geometry[1], cross1 ));
					const b= emit(lerp( geometry[0], geometry[3], cross3 ));
					const c= emit(lerp( geometry[2], geometry[1], cross2 ));
					const d= emit(lerp( geometry[2], geometry[3], cross4 ));

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
                                                 e2([emit(lerp( geometry[2], geometry[1], cross2 )),
		                                         emit(lerp( geometry[0], geometry[1], cross1 )),
		                                         emit(lerp( geometry[3], geometry[1], cross3 ))]);
	                                } else {
                                                 e2([emit(lerp( geometry[0], geometry[1], cross1 )),
	        	                                 emit(lerp( geometry[2], geometry[1], cross2 )),
	                	                         emit(lerp( geometry[3], geometry[1], cross3 ))]);
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
        				const a= emit(lerp( geometry[0], geometry[2], cross1 ));
					const b= emit(lerp( geometry[1], geometry[2], cross2 ));
					const c= emit(lerp( geometry[0], geometry[3], cross3 ));
					const d= emit(lerp( geometry[1], geometry[3], cross4 ));
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
                                                e2( [ emit( lerp( geometry[1], geometry[2], cross2 ) ),
		                                        emit( lerp( geometry[0], geometry[2], cross1 ) ),
		                                        emit( lerp( geometry[3], geometry[2], cross3 ) ) ] );
	                                } else {
                                                e2( [ emit(  lerp( geometry[0], geometry[2], cross1 ) ),
	                                         emit(lerp( geometry[1], geometry[2], cross2 ) ),
	                                         emit(lerp( geometry[3], geometry[2], cross3 ) )]);
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
                                                e2( [ emit( lerp( geometry[0], geometry[3], cross1 ) ),
	                                              emit(  lerp( geometry[1], geometry[3], cross2 ) ),
	                                              emit(  lerp( geometry[2], geometry[3], cross3 ) )]);
	                                } else {
                                                e2( [ emit( lerp( geometry[1], geometry[3], cross2 ) ),
	                                              emit(  lerp( geometry[0], geometry[3], cross1 ) ),
	                                              emit(  lerp( geometry[2], geometry[3], cross3 ) )]);
	                                }
                                }
                        }
                }
        } 
}

// values input to this are in 2 planes for lower and upper values

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



const tetTable = [
	[ 0,2,4,1 ],
	[ 3,1,7,2 ],
	[ 5,7,1,4 ],
	[ 6,7,4,2 ],
	[ 4,7,1,2 ],
]

// 0-5 5-3 
// 3-6 6-0

const tetTable2 = [
	[ 1,0,5,3 ],
	[ 2,3,6,0 ],
	[ 4,5,0,6 ],
	[ 7,3,5,6 ],
	[ 0,6,5,3 ],
]

function cellCompute( alt, values, geom ) {

	_debug && showValues( values );
	_debug && console.log( "Our local state:", values, tets );

	let v = [0,0,0,0];
	let g = [null,null,null,null];
	let n = 0;
	if( alt ) {
		for( let tet of tetTable ) {
			v[0] = values[tet[0]>>2][tet[0]%4];
			v[1] = values[tet[1]>>2][tet[1]%4];
			v[2] = values[tet[2]>>2][tet[2]%4];
			v[3] = values[tet[3]>>2][tet[3]%4];
			g[0] = geom[tet[0]];
			g[1] = geom[tet[1]];
			g[2] = geom[tet[2]];
			g[3] = geom[tet[3]];
	        	tetCompute( v, g, false );
		}
	} else {
		for( let tet of tetTable2 ) {
			v[0] = values[tet[0]>>2][tet[0]%4];
			v[1] = values[tet[1]>>2][tet[1]%4];
			v[2] = values[tet[2]>>2][tet[2]%4];
			v[3] = values[tet[3]>>2][tet[3]%4];
			g[0] = geom[tet[0]];
			g[1] = geom[tet[1]];
			g[2] = geom[tet[2]];
			g[3] = geom[tet[3]];
	        	tetCompute( v, g, false );
		}
	}
}


	const values = [[-1,-1,-1,-1],[-1,-1,-1,-1]];
	for( let x = -1; x < dims[0]; x++ ) {
		cellOrigin[0] = x;
		for( let y = -1; y < dims[1]; y++ ) {
			cellOrigin[1] = y;
			for( let z = 0; z < dims[2]; z++ ) {
				cellOrigin[2] = z;
				const tmp = values[0]; // swap rows
				values[0] = values[1];
				values[1] = tmp;
				
				values[1][0] = ( x < 0 )         ?-1: (0-data[ z * dims[0]*dims[1] + (y+0) * dims[0] + (x+0) ]);
				values[1][1] = ( x >= dims[0]-1 )?-1: (0-data[ z * dims[0]*dims[1] + (y+0) * dims[0] + (x+1) ]);

				if(y  < dims[1]-1 ) {
					values[1][2] = ( x < 0 )?-1: (0-data[ z * dims[0]*dims[1] + (y+1) * dims[0] + (x+0) ]);
					values[1][3] = ( x >= dims[0]-1 )?-1:( 0-data[ z * dims[0]*dims[1] + (y+1) * dims[0] + (x+1) ] );
				} else {
					values[1][2] = -1;
					values[1][3] = -1;
				}
				cellCompute( (x+y+z)&1, values, geom );
			}
		}
	}



  
  return { vertices: vertices, faces: faces };
}
})();


if("undefined" != typeof exports) {
  exports.mesher = MarchingTetrahedra;
}
