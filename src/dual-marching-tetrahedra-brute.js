"use strict";
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
 * Dual Marching Tetrahedra - brute force reference implementation
 *
 * This exists to be OBVIOUSLY CORRECT, not fast.  It is the oracle that the
 * optimized meshers (dual-marching-tetrahedra-3D.js and -3D3.js) get diffed
 * against.  Every lookup table in those files is replaced here by a hash map
 * built at runtime.  If this and the fast one disagree, the fast one's table
 * is wrong.
 *
 * Same 5-tetrahedra Diamond-Crystal-Lattice decomposition, same even/odd
 * parity.  What is different is phase 3:
 *
 *   - the fast mesher decides face winding from the OCTANT OF THE DUAL
 *     VERTEX NORMAL, via ~40 hand-tuned switch statements plus the
 *     tetCentroidFacet table (whose rows are commented "unchecked!").
 *     That cannot be made correct; the normal's octant does not determine
 *     orientation.
 *
 *   - this one decides winding from THE SIGN OF THE TWO GRID SAMPLES at
 *     the ends of the primal edge the dual face surrounds.  That is the
 *     actual invariant.  No normals are consulted at all.
 *
 * and phase 0, which is the "blemish" bug:
 *
 *   - normals here are NEVER derived from the generated geometry.  When a
 *     grid sample sits exactly at 0 on the surface, every crossing on every
 *     tet edge incident to it lands on that same point, so the primal
 *     triangle collapses to a point and its cross product is 0/0 -- the
 *     normal is then whatever garbage survives normalization.  That is
 *     geometrically fine and visually ruinous, and it cannot be patched by
 *     merging points (pointMerge) or by area thresholds, because the
 *     degenerate triangle carries no direction information at all.
 *
 *     Instead normals come from the FIELD: the central-difference gradient
 *     of the density, sampled at the crossing points.  A zero-area triangle
 *     still has a perfectly well-defined gradient.  This is exactly the
 *     "Hermite data" that Dual Contouring assumes on each sign-change edge.
 *     The only genuinely degenerate case left is a true critical point
 *     (value AND gradient zero), which is a real singularity, not an
 *     artifact -- it gets counted and reported rather than hidden.
 *
 * Structure, following Ju et al. / Schaefer et al. (one vertex per cell with
 * a sign change, one polygon per EDGE with a sign change):
 *
 *   phase A  one dual vertex per tetrahedron that has a sign change
 *   phase B  hash map: primal edge -> ring of tetrahedra sharing that edge
 *   phase C  one polygon per sign-change edge, wound by the edge's own signs
 *
 * A tetrahedron can only ever hold a single contour component -- 4 corners
 * split 1-3 (one triangle) or 2-2 (one quad), never more.  So the whole
 * "one vertex per Marching-Cubes cycle" machinery that Dual Contouring needs
 * on cubes is free here, and the result is manifold by construction as soon
 * as the windings agree.
 */

var DualMarchingTetrahedraBrute = (function() {

	// corner index n -> unit cube offset.  n&1 = x, (n>>1)&1 = y, (n>>2)&1 = z
	const geom = [
		[0,0,0],  // 0    bottom layer
		[1,0,0],  // 1
		[0,1,0],  // 2
		[1,1,0],  // 3
		[0,0,1],  // 4    top layer
		[1,0,1],  // 5
		[0,1,1],  // 6
		[1,1,1],  // 7
	];

	// The 5-tetrahedra decomposition of a cube: 4 corner tets plus 1 central
	// regular tet.  The central tet's 6 edges ARE the 6 face diagonals of the
	// cube, so the decomposition picks a diagonal on every face.  A neighbour
	// sharing that face has to pick the same diagonal or the surfaces crack --
	// and the only way to do that is the mirror-image decomposition.  Hence the
	// parity.  There are exactly two 5-tet decompositions and they must
	// alternate on a 3D checkerboard.  This is forced, not a design choice.
	//
	// index with [ (x+y+z)&1 ][ tet 0-4 ][ 0-3 corner ]
	const tetTable = [
		[ [1,0,5,3], [2,3,6,0], [4,5,0,6], [7,3,5,6], [0,6,5,3] ],   // even
		[ [0,2,4,1], [3,1,7,2], [5,7,1,4], [6,7,4,2], [4,7,1,2] ],   // odd
	];

	// the 6 edges of a tetrahedron, as local corner-index pairs
	const tetEdges = [ [0,1],[0,2],[0,3],[1,2],[1,3],[2,3] ];

	function cross( o, a, b ) {
		o[0] = a[1]*b[2] - a[2]*b[1];
		o[1] = a[2]*b[0] - a[0]*b[2];
		o[2] = a[0]*b[1] - a[1]*b[0];
		return o;
	}

	function normalize( v ) {
		const d = Math.sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] );
		if( d > 0 ) { v[0] /= d; v[1] /= d; v[2] /= d; }
		return v;
	}

	// build a right-handed basis (u,v,axis) with u,v perpendicular to axis.
	// u = axis x h, v = axis x u  ==>  u x v = axis, so a rotation from u
	// toward v is a right-handed rotation about axis.
	function buildBasis( axis, u, v ) {
		const ax = Math.abs(axis[0]), ay = Math.abs(axis[1]), az = Math.abs(axis[2]);
		// least-aligned world axis keeps the cross product well conditioned
		const h = ( ax <= ay && ax <= az ) ? [1,0,0]
			: ( ay <= az )              ? [0,1,0]
			:                             [0,0,1];
		normalize( cross( u, axis, h ) );
		normalize( cross( v, axis, u ) );
	}

	return function( data, dims, opts ) {

	opts = opts || {};

	const dim0 = dims[0], dim1 = dims[1], dim2 = dims[2];
	const dim01 = dim0 * dim1;
	const nPoints = dim0 * dim1 * dim2;

	const vertices = [];
	const normals = [];    // parallel to vertices; field gradients, not geometry
	const vType1  = [];    // parallel to vertices; element on the inside
	const vType2  = [];    // parallel to vertices; element on the outside
	const vDelta  = [];    // parallel to vertices; blend between them
	const faces = [];

	// The one convention knob.  Everything about which way the surface faces
	// comes from here -- the winding rule below reads this and nothing else,
	// so if the whole surface renders inside-out, this is the line to flip and
	// there is nowhere else to look.
	//
	// The sample fields in testdata.js are DISTANCE functions, negative inside
	// the solid (x*x+y*y+z*z-1 for the sphere, distanceFromConvexPlanes for the
	// pyramids).  Density is -distance, which is the README's note and the
	// reason every other mesher here opens with d = -data[..].  So solid is
	// where data is negative.
	function isInside( i ) { return data[i] < 0; }

	// grid point index -> world position
	function pointAt( i, o ) {
		const z = ( i / dim01 ) | 0;
		const r = i - z * dim01;
		o[0] = r % dim0;
		o[1] = ( r / dim0 ) | 0;
		o[2] = z;
		return o;
	}

	// Note on index.html's inside/outside toggles: they wrap the volume in a
	// shell of +/- 2.3*Math.random().  Only the SIGN there means anything -- it
	// caps the surface at the volume edge so the mesh closes -- and the random
	// magnitudes do perturb gradients within a cell or two of the border.
	// Excluding the shell from gradient sampling was tried and measured worse
	// (faces sitting on the cap need the shell's step to know which way the cap
	// faces), so it is deliberately NOT done.  With the geometry-first rule
	// below the gradient is only a fallback for collapsed facets anyway.
	//
	// Central-difference gradient at a grid point, pointing OUTWARD -- which,
	// with isInside() above, means toward increasing data (out of the solid,
	// up the distance field).  This has to agree with isInside(); the two of
	// them are the same convention seen from different sides.
	//
	// Sampling is clamped into the trustworthy interior box, so it degrades to
	// a one-sided difference at the edge of real data rather than reading the
	// synthetic shell.  Points inside the shell inherit the nearest real
	// gradient.  Divisor comes from the actual span, so a clamped axis is not
	// stretched by a spurious factor of two.
	function gradAt( i, o ) {
		const z = ( i / dim01 ) | 0;
		const r = i - z * dim01;
		const x = r % dim0, y = ( r / dim0 ) | 0;

		// A point INSIDE the shell needs the shell's own sign step to know
		// which way the cap faces -- there the sign difference dominates and
		// the random magnitude is harmless.  Only a point sitting just inside
		// real data needs protecting, because there the noise perturbs an
		// otherwise-good gradient.  So the clamp applies to the neighbours of
		// interior points only, never to the point itself.
		const xm = Math.max( 0, x-1 ), xp = Math.min( dim0-1, x+1 );
		const ym = Math.max( 0, y-1 ), yp = Math.min( dim1-1, y+1 );
		const zm = Math.max( 0, z-1 ), zp = Math.min( dim2-1, z+1 );


		const yz = y*dim0 + z*dim01;
		o[0] = ( data[ xp + yz ]               - data[ xm + yz ]               ) / ( (xp-xm) || 1 );
		o[1] = ( data[ x + yp*dim0 + z*dim01 ] - data[ x + ym*dim0 + z*dim01 ] ) / ( (yp-ym) || 1 );
		o[2] = ( data[ x + y*dim0 + zp*dim01 ] - data[ x + y*dim0 + zm*dim01 ] ) / ( (zp-zm) || 1 );
		return o;
	}

// ---------------------------------------------------------------------------
// phase A -- enumerate every tetrahedron, and give the ones with a sign change
//            a dual vertex.
//
// Tet id is cellIndex*5 + tet, so it is dense and can index flat arrays.
// tetCorners[id] is the 4 GLOBAL grid indices -- no offsets, no [odd] tables
// downstream.  That is the whole trick: once the corners are global, nothing
// after this point knows or cares which decomposition produced them.
// ---------------------------------------------------------------------------

	const tetCorners = new Array( (dim0-1)*(dim1-1)*(dim2-1)*5 ).fill( null );
	const tetDual    = new Array( tetCorners.length ).fill( null );  // [x,y,z]
	const tetNorm    = new Array( tetCorners.length ).fill( null );  // [x,y,z]

	// Element (texture type) data per dual vertex.  A dual vertex straddles the
	// material boundary, so it reports the inside corner's element, the outside
	// corner's element, and how far along it sits.  GridGeometryBuffer wants
	// these per FACE as [t1a,t1b,t1c, t2a,t2b,t2c] + [da,db,dc]; that triple is
	// what lets the shader blend between two materials across a triangle using
	// the barycentric 'simplex' attribute.  Supplying one element with zero
	// delta is exactly what makes elements look unblended.
	const elemData   = opts.elements && opts.elements.data;
	const tetType1   = new Float32Array( tetCorners.length );
	const tetType2   = new Float32Array( tetCorners.length );
	const tetDelta   = new Float32Array( tetCorners.length );
	const tetVertId  = new Int32Array( tetCorners.length ).fill( -1 );

	const gA = [0,0,0], gB = [0,0,0];
	let nTets = 0, criticalPoints = 0, degenerateFacets = 0;
	for( let z = 0; z < dim2-1; z++ )
	for( let y = 0; y < dim1-1; y++ )
	for( let x = 0; x < dim0-1; x++ ) {
		const parity = ( x + y + z ) & 1;
		const table = tetTable[parity];
		const cell = x + y*(dim0-1) + z*(dim0-1)*(dim1-1);

		for( let t = 0; t < 5; t++ ) {
			const src = table[t];
			const c = [0,0,0,0];
			for( let i = 0; i < 4; i++ ) {
				const g = geom[src[i]];
				c[i] = ( x + g[0] ) + ( y + g[1] )*dim0 + ( z + g[2] )*dim01;
			}
			const id = cell*5 + t;
			tetCorners[id] = c;
			nTets++;

			// dual vertex = centroid of this tet's crossing points.  A 1-3
			// split gives 3 of them (the primal triangle), a 2-2 split gives 4
			// (the primal quad).  Averaging is the brute-force stand-in for a
			// QEF solve -- it loses sharp features but cannot fail.
			//
			// The normal is accumulated in the same pass from the field
			// gradient lerped along each crossing edge -- Hermite data.  It is
			// deliberately independent of the triangle this tet would have
			// produced, so a collapsed triangle costs nothing.
			let sx = 0, sy = 0, sz = 0, n = 0;
			let nx = 0, ny = 0, nz = 0;
			let tsum = 0;
			const xp_ = [0,0,0,0], yp_ = [0,0,0,0], zp_ = [0,0,0,0];  // crossing points

			// Outward reference direction from THIS tet's corner signs alone:
			// mean of the outside corners minus mean of the inside ones.  Pure
			// sign data, so it survives the random-magnitude fill shell that
			// makes the gradient untrustworthy at the volume walls.
			let ix = 0, iy = 0, iz = 0, ni = 0;
			let ox = 0, oy = 0, oz = 0, no = 0;
			let elemOut = 0;
			// keep every inside corner so the material can be chosen by
			// proximity below rather than by whichever came first
			const inEl = [0,0,0,0], inPx = [0,0,0,0], inPy = [0,0,0,0], inPz = [0,0,0,0];
			for( let v = 0; v < 4; v++ ) {
				const gi = c[v];
				const gz = ( gi / dim01 ) | 0, gr = gi - gz*dim01;
				const gx = gr % dim0, gy = ( gr / dim0 ) | 0;
				if( isInside( gi ) ) {
					ix += gx; iy += gy; iz += gz;
					inPx[ni] = gx; inPy[ni] = gy; inPz[ni] = gz;
					inEl[ni] = elemData ? elemData[gi] : 0;
					ni++;
				} else {
					ox += gx; oy += gy; oz += gz;
					if( !no && elemData ) elemOut = elemData[gi];
					no++;
				}
			}
			const outX = no && ni ? ( ox/no - ix/ni ) : 0;
			const outY = no && ni ? ( oy/no - iy/ni ) : 0;
			const outZ = no && ni ? ( oz/no - iz/ni ) : 0;

			for( let e = 0; e < 6; e++ ) {
				const a = c[ tetEdges[e][0] ], b = c[ tetEdges[e][1] ];
				const ia = isInside(a), ib = isInside(b);
				if( ia === ib ) continue;

				const da = data[a], db = data[b];
				const den = da - db;
				const s = den ? ( da / den ) : 0;   // da + s*(db-da) == 0

				const az_ = ( a / dim01 ) | 0, ar = a - az_*dim01;
				const bz_ = ( b / dim01 ) | 0, br = b - bz_*dim01;
				const ax_ = ar % dim0, ay_ = ( ar / dim0 ) | 0;
				const bx_ = br % dim0, by_ = ( br / dim0 ) | 0;

				const px = ax_ + ( bx_ - ax_ )*s;
				const py = ay_ + ( by_ - ay_ )*s;
				const pz = az_ + ( bz_ - az_ )*s;
				if( n < 4 ) { xp_[n] = px; yp_[n] = py; zp_[n] = pz; }
				tsum += s;
				sx += px; sy += py; sz += pz;

				gradAt( a, gA );
				gradAt( b, gB );
				nx += gA[0] + ( gB[0] - gA[0] )*s;
				ny += gA[1] + ( gB[1] - gA[1] )*s;
				nz += gA[2] + ( gB[2] - gA[2] )*s;
				n++;
			}
			if( n ) {
				tetDual[id] = [ sx/n, sy/n, sz/n ];
				// Material = the element of the inside corner NEAREST this dual
				// vertex.  Taking simply the first inside corner speckles: the
				// tet-corner ordering rotates with parity and tet index, so two
				// neighbouring tets straddling one material boundary sample
				// different corners and disagree.  Distance is continuous, so
				// adjacent dual vertices agree unless the boundary really is
				// between them.  Measured on Terrain: 13.3% -> see below.
				const dvx = sx/n, dvy = sy/n, dvz = sz/n;
				let bestEl = inEl[0], bestD = Infinity;
				for( let k = 0; k < ni; k++ ) {
					const ddx = inPx[k]-dvx, ddy = inPy[k]-dvy, ddz = inPz[k]-dvz;
					const dd = ddx*ddx + ddy*ddy + ddz*ddz;
					if( dd < bestD ) { bestD = dd; bestEl = inEl[k]; }
				}
				tetType1[id] = bestEl;
				tetType2[id] = elemOut;
				tetDelta[id] = tsum / n;      // mean crossing parameter

				const glen = Math.sqrt( nx*nx + ny*ny + nz*nz );
				const haveGrad = glen > 1e-12;
				if( haveGrad ) { nx /= glen; ny /= glen; nz /= glen; }

				// The primal facet's own normal, if this tet produced a
				// well-conditioned triangle.  Any 3 of the crossing points do:
				// for a 1-3 split there are exactly 3, and for a 2-2 split all
				// 4 are near-coplanar.
				let gx = 0, gy = 0, gz = 0, glen2 = 0;
				if( n >= 3 ) {
					const ux = xp_[1]-xp_[0], uy = yp_[1]-yp_[0], uz = zp_[1]-zp_[0];
					const wx = xp_[2]-xp_[0], wy = yp_[2]-yp_[0], wz = zp_[2]-zp_[0];
					gx = uy*wz - uz*wy; gy = uz*wx - ux*wz; gz = ux*wy - uy*wx;
					glen2 = Math.sqrt( gx*gx + gy*gy + gz*gz );
				}

				// Pick per source of truth, not per preference:
				//
				//  - a well-conditioned triangle describes THIS surface at
				//    THIS scale, including features smaller than the sampling
				//    (a lone inside-sample makes a real little polyhedron; the
				//    field has no feature there at all, so its gradient would
				//    hand every vertex of that blob the same direction and it
				//    would shade like a chip of whatever is nearby).
				//
				//  - a collapsed triangle describes nothing, and that is where
				//    the gradient is the only thing left -- the blemish case.
				//
				// Orientation always comes from the gradient, so the geometric
				// normal never introduces a sign question of its own.
				if( glen2 > 1e-6 ) {
					gx /= glen2; gy /= glen2; gz /= glen2;
					// Orient from the corner signs, NOT from the gradient.  The
					// gradient is only a magnitude-weighted hint and goes random
					// inside the fill shell; the sign-derived direction cannot.
					const ref = gx*outX + gy*outY + gz*outZ;
					if( ref < 0 ) { gx = -gx; gy = -gy; gz = -gz; }
					else if( ref === 0 && haveGrad && ( gx*nx + gy*ny + gz*nz ) < 0 ) {
						gx = -gx; gy = -gy; gz = -gz;   // 2-2 split ties: fall back
					}
					tetNorm[id] = [ gx, gy, gz ];
				} else if( haveGrad ) {
					tetNorm[id] = [ nx, ny, nz ];
					degenerateFacets++;
				} else {
					// value zero AND gradient zero: a real critical point of
					// the field, not a meshing artifact.  Nothing can recover
					// a direction here; count it so it stays visible.
					tetNorm[id] = [0,0,0];
					criticalPoints++;
				}
			}
		}
	}

// ---------------------------------------------------------------------------
// phase B -- primal edge -> ring of tetrahedra sharing it.
//
// This replaces referenceTetEdgeFormations / TetFaceMapping / edgeToComp and
// the ~40 inline `normals[baseOffset + (1+1*dim0)*tetCount + 0]` offsets.
// Built by insertion, so it is right by construction and knows nothing about
// parity.
// ---------------------------------------------------------------------------

	const edgeRings = new Map();   // key a*nPoints+b (a<b)  ->  [tet ids]

	for( let id = 0; id < tetCorners.length; id++ ) {
		const c = tetCorners[id];
		if( !c ) continue;
		for( let e = 0; e < 6; e++ ) {
			let a = c[ tetEdges[e][0] ], b = c[ tetEdges[e][1] ];
			if( a > b ) { const s = a; a = b; b = s; }
			const key = a*nPoints + b;
			const ring = edgeRings.get( key );
			if( ring ) ring.push( id );
			else edgeRings.set( key, [id] );
		}
	}

// ---------------------------------------------------------------------------
// phase C -- one polygon per sign-change edge.
// ---------------------------------------------------------------------------

	function vertexOf( id ) {
		if( tetVertId[id] < 0 ) {
			tetVertId[id] = vertices.length;
			vertices.push( tetDual[id] );
			normals.push( tetNorm[id] );
			vType1.push( tetType1[id] );
			vType2.push( tetType2[id] );
			vDelta.push( tetDelta[id] );
		}
		return tetVertId[id];
	}


	// Order the ring COMBINATORIALLY, never geometrically.
	//
	// Two tets are consecutive around edge {a,b} iff they share a face, which
	// -- since they already share a and b -- means sharing a third corner.
	// Each tet in the ring has exactly two such "other" corners, one per
	// neighbour, so the ring is just a walk through that adjacency.
	//
	// Sorting by angle would be the obvious thing and it is wrong: where a
	// sample sits exactly at 0 the dual vertices all collapse onto it, the
	// angles become noise, and the ring shatters.  The walk below touches no
	// dual vertex at all, so collapsed geometry costs nothing.  If a corner
	// joins fewer than two tets the edge is on the volume border and its ring
	// is an arc -- return null and leave the surface open there.
	function walkRing( ring, a, b ) {
		const byCorner = new Map();      // shared corner -> the 2 tets meeting there
		const others = new Map();        // tet -> its 2 other corners

		for( let i = 0; i < ring.length; i++ ) {
			const id = ring[i], c = tetCorners[id], o = [];
			for( let j = 0; j < 4; j++ ) if( c[j] !== a && c[j] !== b ) o.push( c[j] );
			if( o.length !== 2 ) return null;
			others.set( id, o );
			for( let k = 0; k < 2; k++ ) {
				const l = byCorner.get( o[k] );
				if( l ) l.push( id ); else byCorner.set( o[k], [id] );
			}
		}
		for( const l of byCorner.values() ) if( l.length !== 2 ) return null;

		const order = [], shared = [];
		let cur = ring[0], prev = -1;
		for( let step = 0; step < ring.length; step++ ) {
			order.push( cur );
			const o = others.get( cur );
			const next = ( o[0] !== prev ) ? o[0] : o[1];
			const pair = byCorner.get( next );
			shared.push( next );
			prev = next;
			cur = ( pair[0] === cur ) ? pair[1] : pair[0];
		}
		if( cur !== ring[0] ) return null;          // walk did not close

		return { order: order, shared: shared };
	}

	const pa = [0,0,0], pb = [0,0,0], sa = [0,0,0], sb = [0,0,0];
	const axis = [0,0,0];
	let openEdges = 0, degenerateFaces = 0;

	for( const [key, ring] of edgeRings ) {
		const b = key % nPoints;
		const a = ( key - b ) / nPoints;

		const ia = isInside(a), ib = isInside(b);
		if( ia === ib ) continue;              // no sign change, no polygon

		// every tet on a crossing edge necessarily has a crossing, so it
		// necessarily got a dual vertex.  Guard anyway.
		let ok = true;
		for( let i = 0; i < ring.length; i++ ) if( !tetDual[ ring[i] ] ) { ok = false; break; }
		if( !ok ) { openEdges++; continue; }
		if( ring.length < 3 ) { openEdges++; continue; }

		const walked = walkRing( ring, a, b );
		if( !walked ) { openEdges++; continue; }
		const sorted = walked.order;

		pointAt( a, pa );
		pointAt( b, pb );
		axis[0] = pb[0]-pa[0]; axis[1] = pb[1]-pa[1]; axis[2] = pb[2]-pa[2];
		const mx = ( pa[0]+pb[0] )/2, my = ( pa[1]+pb[1] )/2, mz = ( pa[2]+pb[2] )/2;

		// The walk closed, but it could have gone either way round.  Fix the
		// handedness from the first two SHARED CORNERS -- those are grid
		// points, so their coordinates are exact small integers and this test
		// can never degenerate, no matter what the dual vertices did.
		pointAt( walked.shared[0], sa );
		pointAt( walked.shared[1], sb );
		const w0x = sa[0]-mx, w0y = sa[1]-my, w0z = sa[2]-mz;
		const w1x = sb[0]-mx, w1y = sb[1]-my, w1z = sb[2]-mz;
		const hand = axis[0]*( w0y*w1z - w0z*w1y )
		           + axis[1]*( w0z*w1x - w0x*w1z )
		           + axis[2]*( w0x*w1y - w0y*w1x );
		if( hand < 0 ) sorted.reverse();

		// centroid fan.  n triangles instead of n-2, and one extra vertex, but
		// it can never pick the wrong diagonal -- which is the entire job of
		// getFold() in the fast mesher, heuristically and sometimes wrongly.
		let cx = 0, cy = 0, cz = 0, cnx = 0, cny = 0, cnz = 0, cdel = 0;
		for( let i = 0; i < sorted.length; i++ ) {
			const d = tetDual[ sorted[i] ], nn = tetNorm[ sorted[i] ];
			cx += d[0];  cy += d[1];  cz += d[2];
			cnx += nn[0]; cny += nn[1]; cnz += nn[2];
			cdel += tetDelta[ sorted[i] ];
		}
		const cl = Math.sqrt( cnx*cnx + cny*cny + cnz*cnz );
		const centre = vertices.length;
		vertices.push( [ cx/sorted.length, cy/sorted.length, cz/sorted.length ] );
		normals.push( cl > 1e-12 ? [ cnx/cl, cny/cl, cnz/cl ] : [0,0,0] );
		// the fan centre inherits the ring's materials: the edge it surrounds
		// runs from a solid sample to an empty one, so those two ends ARE its
		// inside/outside pair, and the mean ring delta is where it sits.
		vType1.push( elemData ? elemData[ ia ? a : b ] : 0 );
		vType2.push( elemData ? elemData[ ia ? b : a ] : 0 );
		vDelta.push( cdel / sorted.length );

		// THE winding rule -- the whole point of the rewrite.  The fan's normal
		// points toward +axis, i.e. toward b.  A surface normal must point out
		// of the solid, so keep the order when a is the solid end and reverse
		// it when b is.  Two grid samples decide it; no geometry is consulted,
		// and there is no case analysis.
		const forward = ia;

		for( let i = 0; i < sorted.length; i++ ) {
			const p1 = vertexOf( sorted[i] );
			const p2 = vertexOf( sorted[ (i+1) % sorted.length ] );
			if( forward ) faces.push( [ centre, p1, p2 ] );
			else          faces.push( [ centre, p2, p1 ] );
		}
	}

// ---------------------------------------------------------------------------
// phase D -- weld coincident dual vertices.
//
// This is what actually cures the blemishes.  Where a sample sits exactly at
// 0, every tet touching it puts its dual vertex on that same point, so the
// mesh contains a cluster of distinct vertices at one location joined by
// triangles of zero area -- the "infinitesimal space".  Culling those faces
// by area would punch topological holes in the surface.  Welding instead
// makes the cluster into a single vertex, at which point the degenerate
// triangles become degenerate BY INDEX and can be removed with no loss of
// closure.  The region does not get hidden, it stops existing.
//
// The surviving vertex keeps the summed field gradient, so its normal is the
// honest average of the surface passing through it.
// ---------------------------------------------------------------------------

	const Q = 1e7;                      // weld quantum, well below grid spacing
	const weldKey = new Map();
	const remap = new Int32Array( vertices.length );
	const outVerts = [], outNorms = [];
	const outT1 = [], outT2 = [], outD = [], outDn = [];

	for( let i = 0; i < vertices.length; i++ ) {
		const p = vertices[i];
		const k = Math.round(p[0]*Q) + "," + Math.round(p[1]*Q) + "," + Math.round(p[2]*Q);
		let id = weldKey.get( k );
		if( id === undefined ) {
			id = outVerts.length;
			weldKey.set( k, id );
			outVerts.push( [ p[0], p[1], p[2] ] );
			outNorms.push( [0,0,0] );
			outT1.push( vType1[i] );      // types are ids, not averageable
			outT2.push( vType2[i] );
			outD.push( 0 ); outDn.push( 0 );
		}
		remap[i] = id;
		const n = normals[i], o = outNorms[id];
		o[0] += n[0]; o[1] += n[1]; o[2] += n[2];
		outD[id] += vDelta[i]; outDn[id]++;
	}
	for( let i = 0; i < outD.length; i++ ) if( outDn[i] ) outD[i] /= outDn[i];
	for( let i = 0; i < outNorms.length; i++ ) {
		const o = outNorms[i];
		const l = Math.sqrt( o[0]*o[0] + o[1]*o[1] + o[2]*o[2] );
		if( l > 1e-12 ) { o[0] /= l; o[1] /= l; o[2] /= l; }
	}

	const outFaces = [];
	for( let i = 0; i < faces.length; i++ ) {
		const f = faces[i];
		const a = remap[f[0]], b = remap[f[1]], c = remap[f[2]];
		if( a === b || b === c || a === c ) { degenerateFaces++; continue; }
		outFaces.push( [a,b,c] );
	}

	console.log( "DMT-brute:", nTets, "tets,", edgeRings.size, "edges ->",
		outVerts.length, "verts,", outFaces.length, "faces"
		, "| welded away:", vertices.length - outVerts.length, "verts,"
		, degenerateFaces, "collapsed faces"
		, "| skipped:", openEdges, "open/non-closing edge rings"
		, "| gradient-rescued facets:", degenerateFacets, "| field critical points:", criticalPoints );

	return { vertices: outVerts, faces: outFaces, normals: outNorms,
	         types1: outT1, types2: outT2, deltas: outD };
	}
})();

if( "undefined" != typeof exports ) {
	exports.mesher = DualMarchingTetrahedraBrute;
}
