// Seedable random number generators, so test data and lights are the same
// every run (Math.random() made every screenshot different).
//
// Implemented from: https://stackoverflow.com/questions/521295/seeding-the-random-number-generator-in-javascript
// more https://github.com/bryc/code/blob/master/jshash/PRNGs.md#lcg-lehmer-rng
//
//   const rand = SFC32( "apples" );  // seed from any string
//   rand();                          // 0 <= n < 1, like Math.random()

// used to get a seed from a string.
function xmur3(str) {
    for(var i = 0, h = 1779033703 ^ str.length; i < str.length; i++) {
        h = Math.imul(h ^ str.charCodeAt(i), 3432918353);
        h = h << 13 | h >>> 19;
    } return function() {
        h = Math.imul(h ^ (h >>> 16), 2246822507);
        h = Math.imul(h ^ (h >>> 13), 3266489909);
        return (h ^= h >>> 16) >>> 0;
    }
}

// simple fast counter
// sfc32 is part of the PractRand random number testing suite (which it passes of course).
// sfc32 has a 128-bit state and is very fast in JS.
function sfc32(a, b, c, d) {
    return function() {
      a >>>= 0; b >>>= 0; c >>>= 0; d >>>= 0; 
      var t = (a + b) | 0;
      a = b ^ b >>> 9;
      b = c + (c << 3) | 0;
      c = (c << 21 | c >>> 11);
      d = d + 1 | 0;
      t = t + d | 0;
      c = c + t | 0;
      return (t >>> 0) / 4294967296;
    }
}

// Mulberry32 is a simple generator with a 32-bit state, but is extremely fast and has
// good quality randomness (author states it passes all tests of gjrand testing suite
// and has a full 2^32 period).
function mulberry32(a) {
    return function() {
      var t = a += 0x6D2B79F5;
      t = Math.imul(t ^ t >>> 15, t | 1);
      t ^= t + Math.imul(t ^ t >>> 7, t | 61);
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    }
}

// xoshiro128** (Vigna & Blackman, 2018), the fastest generator with a 128-bit state.
function xoshiro128ss(a, b, c, d) {
    return function() {
        var t = b << 9, r = a * 5; r = (r << 7 | r >>> 25) * 9;
        c ^= a; d ^= b;
        b ^= c; a ^= d; c ^= t;
        d = d << 11 | d >>> 21;
        return (r >>> 0) / 4294967296;
    }
}

// JSF or 'smallprng' by Bob Jenkins (2007).  Passes PractRand; not as fast as sfc32.
function jsf32(a, b, c, d) {
    return function() {
        a |= 0; b |= 0; c |= 0; d |= 0;
        var t = a - (b << 27 | b >>> 5) | 0;
        a = b ^ (c << 17 | c >>> 15);
        b = c + d | 0;
        c = d + t | 0;
        d = a + t | 0;
        return (d >>> 0) / 4294967296;
    }
}

export function getSeed(s) { return xmur3( s ); }
export function getSeed4(seed) { const s = xmur3( seed ); return [s(),s(),s(),s()]; }

export function SFC32( seed ) { const s = xmur3( seed ); return sfc32(s(),s(),s(),s()); }
export function MUL32( seed ) { const s = xmur3( seed ); return mulberry32(s()); }
export function XOR32( seed ) { const s = xmur3( seed ); return xoshiro128ss(s(),s(),s(),s()); }
export function JSF32( seed ) { const s = xmur3( seed ); return jsf32(s(),s(),s(),s()); }
export function SFC32_( seed ) { return sfc32(...seed); }
export function MUL32_( seed ) { return mulberry32(seed[0]); }
export function XOR32_( seed ) { return xoshiro128ss(...seed); }
export function JSF32_( seed ) { return jsf32(...seed); }
