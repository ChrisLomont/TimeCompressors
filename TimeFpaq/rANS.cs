// based on 
// Simple byte-aligned rANS encoder/decoder - public domain - Fabian 'ryg' Giesen 2014
// 
// C# conversion, lots of rewrite by Chris Lomont, 2022
using System.Diagnostics;

namespace TimeFpaq
{
    internal class rANS
    {

        /*

        need to characterize values in here:
        1. L = 1<<23 lower bound (based on byte 8 bit b output)
        2. prob bits 14, needs to divide L
        3. M = 1<<probbits, M | L
        4. prob scaling, baed on block size when encoding
        5. buffer size - must have limits?
        6. internal word size:L 32 bits, needs 32*32 mult into 64 to be fast
        7. header based on these (todo - store orig freqs, smaller? not scaled? - or scaled smaller?)
        8. reciprocal format and papers
        9. characterize 8 bit (23+8 format) and 16 bit (15+16 format) output



        byte oriented rANS compressor/decompressor - Chris Lomont 2022


        Idea works like this:
        Suppose you have symbols
        TODO - write up

        TODO - integrate oder 0 or 1 model https://github.com/jkbonfield/rans_static
        TODO - https://github.com/jkbonfield/rans_static talks about 15+16 size model instead of this 23+8, allows one renormalization, never 2
         
        TODO - replace Alverson reciprocal division with with Robison 2005 - "N-Bit Unsigned Division Via N-Bit Multiply-Add"
        or maybe Drane 2012 "Correctly Rounded Constant Integer Division via Multiply-Add"

        TODO - incorporate https://github.com/loxxous/rans-tricks
        https://encode.su/threads/2944-Rans-tricks

        geisen posted faster version at https://gist.github.com/rygorous/b434cf2916be5c9573796b5f96671cbe

        header size optimization https://encode.su/threads/1883-Reducing-the-header-optimal-quantization-compression-of-probability-distribution

        fast notes
        http://cbloomrants.blogspot.com/2014/02/02-18-14-ansfast-implementation-notes.html

        */


        // Simple byte-aligned rANS encoder/decoder - public domain - Fabian 'ryg' Giesen 2014
        //
        // Not intended to be "industrial strength"; just meant to illustrate the general
        // idea.

        // Chris Lomont port to C#, rewrite and nicer API for C#, 2022


        // READ ME FIRST:
        //
        // This is designed like a typical arithmetic coder API, but there's three
        // twists you absolutely should be aware of before you start hacking:
        //
        // 1. You need to encode data in *reverse* - last symbol first. rANS works
        //    like a stack: last in, first out.
        // 2. Likewise, the encoder outputs bytes *in reverse* - that is, you give
        //    it a pointer to the *end* of your buffer (exclusive), and it will
        //    slowly move towards the beginning as more bytes are emitted.
        // 3. Unlike basically any other entropy coder implementation you might
        //    have used, you can interleave data from multiple independent rANS
        //    encoders into the same bytestream without any extra signaling;
        //    you can also just write some bytes by yourself in the middle if
        //    you want to. This is in addition to the usual arithmetic encoder
        //    property of being able to switch models on the fly. Writing raw
        //    bytes can be useful when you have some data that you know is
        //    incompressible, and is cheaper than going through the rANS encode
        //    function. Using multiple rANS coders on the same byte stream wastes
        //    a few bytes compared to using just one, but execution of two
        //    independent encoders can happen in parallel on superscalar and
        //    Out-of-Order CPUs, so this can be *much* faster in tight decoding
        //    loops.
        //
        //    This is why all the rANS functions take the write pointer as an
        //    argument instead of just storing it in some context struct.

        // --------------------------------------------------------------------------

        // L ('l' in the paper) is the lower bound of our normalization interval.
        // Between this and our byte-aligned emission, we use 31 (not 32!) bits.
        // This is done intentionally because exact reciprocals for 31-bit uints
        // fit in 32-bit uints: this permits some optimizations during encoding.

        /* Chris Lomont 2022
        Asymmetric Numerical System (ANS) Compression 
        rANS flavor

        Theory:
        based on https://fgiesen.wordpress.com/2014/02/02/rans-notes/
        rewritten as I understand it

        Want to compress symbols in alphabet {s_0, s_1, .., s_{n-1}}.
        Frequency count of symbol i is f_i, F = sum_i f_i is total count
        thus prob of symbol s_i is p_i = f_i/F. 
        Let t_i = sum_j=0^{j=i-1} s_j be the tally symbol frequencies up to 

        Let L = the stream as F symbols, sorted, so symbol s_i occurs f_i in a row. Then let        
        g(x) = symbol in x-th position in L.
        
        Thus for symbol counts (2,3,1,0,5) you have sequence (0,0, 1,1,1, 2, 4,4,4,4,4)
        g(6) = 4, g(0)=1, g(2) = 1, etc.

        Consider a compressor C that takes an integer state x and symbol s_i, and encodes as
        C(s_i,x) = x' = F * Floor(x/f_i) + t_i + (x % f_i)
        and a decompressor D decodes as
        D(x') = (s_i,x) = (g(x' % F), f_i Floor(x'/T) - t_i + (x % F))

        E can be viewed as: 
            1) take the state, floor divide by f_i to an integer
            2) multiply by F, making "space" to store an offset allowing extracting s_i and x again
            3) add b_i, the space made by F, where there are now f_i places to store a little more info
            4) add (x % f_i), capturing the information lost by step 1 which took a floor
        D is the reverse process, taking a state x', extracting a symbol s_i, then using knowledge of s_i to reverse E.
        Note E and D are double sided inverses. 

        This provides an entropy encoder with compression ratio like arithmetic encoder, since a symbol
        with probability p_i = f_i/F increases the state size by ~ -log_2(p_i) bits as expected.

        To make streaming, we need to output bits from x when suitable, and rescale it. 
        Define the normalized interval where we want x to be as I = {L, L+1, ..., bL-1} = [L:bL)
        where b is the size of items we output. b = 2 is bitwise input/output, b=256 is bytewise, etc.

        TODO - see link above for more notes.... finish this note

        This leads to encoder C(x,s) (renormalize, then encode):

        
        while (!done) 
            assert(L <= x && x < b*L) # loop invariant
        
            x_max = (b * (L / F)) * freq[s]; 
            while (x >= x_max) 
                writeToStream(x % b)
                x /= b
            x = freq[s] * (x / F) + (x % F) + base[s]

        divisions can be replaced with shift for fixed powers of 2, 
        the 
        
        The last thing is to note that encoding acts like placing digits 
        at the high low end of a number, and decoding is stripping them off, 
        so one works in reverse to the other. Thus we choose encoding is
        backwards from a buffer, decoding is forwards, spitting out symbols
        in the original order.

        Needs: to ensure above works, need
        L = kF for some positive const k
        F = power of 2 for fast division and modulus
        div by f_i replaced with reciprocal multiplication (*)


        (*) Alvarez TODO 31 bits, or TODO for 32 bits

            // compute floor(x/d) for x in [0,2^N-1], fixed d, using form
            // Floor[(a*x+b)/2^s]
            // given d, compute a,b,s
            (uint a, uint b, int s) Div2(uint d)
            { // Robison, 2005 "N-Bit Unsigned Division Via N-Bit Multiply-Add"

var m = (int)Math.Floor(Math.Log2(d));
                uint a, b;
                int n = 32; // 32 bit integers
                if (d == 1 << m)
                {
                    a = b = UInt32.MaxValue;
                }
                else
                {
                    var t = (uint)((1UL << (m + n)) / d); // as 64 bit
                    var r = ((t * d + d));// & (1UL << (n - 1));
                    if (r <= (1UL << m))
                    {
                        a = (uint)t + 1;
                        b = 0;
                    }
                    else
                    {
                        a = b = (uint)t;
                    }

                }
                return (a, b, m + 32);
        */


        const uint LowerBound = 1u << 23; // lower bound L of our normalization interval
        const uint ProbabilityBits = 14;
        const uint ProbabilityScaling = 1 << (int)ProbabilityBits;

        class SymbolStats
        {
            // gather stats on data
            public void ComputeSymbolStatistics(IEnumerable<byte> data)
            {
                // count frequencies of all symbols
                for (var i = 0; i < 256; i++)
                    Frequencies[i] = 0;

                foreach (var t in data)
                    Frequencies[t]++;

                NormalizeFrequencies();
            }
            
            public uint[] Frequencies { get;  }= new uint[256];

            #region Implementation

            /// <summary>
            /// Scales frequencies to make sum of them be a certain power of 2 (ProbabilityScaling)
            /// </summary>
            void NormalizeFrequencies()
            {

                Debug.Assert(ProbabilityScaling >= 256);

                // todo - rederive, want to minimize loss due to rounding errors, analyze

                // todo - rewrite to avoid need to create/store cumulative freqs
                uint[] cumulativeFrequencies = new uint[257];

                // compute cumulative frequencies
                cumulativeFrequencies[0] = 0;
                for (var i = 0; i < 256; i++)
                    cumulativeFrequencies[i + 1] = cumulativeFrequencies[i] + Frequencies[i];

                var currentTotal = cumulativeFrequencies[256];

                // resample distribution based on cumulative freqs
                for (var i = 1; i <= 256; i++)
                    cumulativeFrequencies[i] = (uint)((ulong)ProbabilityScaling * cumulativeFrequencies[i] / currentTotal);

                // if we nuked any non-0 frequency symbol to 0, we need to steal
                // the range to make the frequency nonzero from elsewhere.
                //
                // this is not at all optimal, i'm just doing the first thing that comes to mind.
                for (var i = 0; i < 256; i++)
                {
                    if (Frequencies[i] != 0 && cumulativeFrequencies[i + 1] == cumulativeFrequencies[i])
                    {
                        // Console.WriteLine($"Zero freq {i}");
                        // symbol i was set to zero freq

                        // find best symbol to steal frequency from (try to steal from low-freq ones)
                        uint bestFrequency = uint.MaxValue;
                        int bestSteal = -1;
                        for (var j = 0; j < 256; j++)
                        {
                            var freq = cumulativeFrequencies[j + 1] - cumulativeFrequencies[j];
                            if (freq > 1 && freq < bestFrequency)
                            {
                                bestFrequency = freq;
                                bestSteal = j;
                            }
                        }

                        Debug.Assert(bestSteal != -1);

                        // and steal from it!
                        if (bestSteal < i)
                        {
                            for (var j = bestSteal + 1; j <= i; j++)
                                cumulativeFrequencies[j]--;
                        }
                        else
                        {
                            Debug.Assert(bestSteal > i);
                            for (var j = i + 1; j <= bestSteal; j++)
                                cumulativeFrequencies[j]++;
                        }
                    }
                }

                // calculate updated freqs and make sure we didn't screw anything up
                Debug.Assert(cumulativeFrequencies[0] == 0 && cumulativeFrequencies[256] == ProbabilityScaling);
                for (var i = 0; i < 256; i++)
                {
                    if (Frequencies[i] == 0)
                        Debug.Assert(cumulativeFrequencies[i + 1] == cumulativeFrequencies[i]);
                    else
                        Debug.Assert(cumulativeFrequencies[i + 1] > cumulativeFrequencies[i]);

                    // calc updated freq
                    Frequencies[i] = cumulativeFrequencies[i + 1] - cumulativeFrequencies[i];
                }
            }
            #endregion
        }


        class Compressor
        {
            public Compressor()
            {
                for (var i = 0; i < 256; i++)
                    symbols[i] = new();
            }

            /*
             Block format:
            byte version (implies parameters?) - only front of all blocks
            block size (4 bytes) - only front of all blocks

            Per block:
                byte size of freqs
                freqs (256 of them, little endian), 
                block bytes of data (may be bigger if end overlapped, up to ?? bytes??)

             */

            public List<byte> Compress(IReadOnlyList<byte> data)
            {
                // todo - make block sized compressor

                symbolStats.ComputeSymbolStatistics(data);

                uint cumulativeFreq = 0;
                for (var i = 0; i < 256; i++)
                {
                    InitializeSymbol(symbols[i], cumulativeFreq, symbolStats.Frequencies[i]);
                    cumulativeFreq += symbolStats.Frequencies[i];
                }

                var compressedData = new List<byte>();

                // compressor state - here to allow multi-stream encoding/decoding
                // state is only a 32 bit integer
                uint state = LowerBound;

                // process data in reverse to compress it
                for (var i = data.Count; i > 0; i--)
                    CompressSymbol(ref state, compressedData, symbols[data[i - 1]]);
                FlushCompressor(state, compressedData);


                // write header, backwards, at end before reverse:
                // todo - investigate shortening this, packing, etc...
                for (var i = 256; i > 0; --i)
                    Write4(symbolStats.Frequencies[i - 1], compressedData);
                compressedData.Add(4); // 1 byte size of freq entry
                Write4((uint)data.Count, compressedData);
                compressedData.Add(0); // version

                // can do above to write backwards if fixed size arrays, then export them as blocks to somewhere
                compressedData.Reverse();
                return compressedData; // costly, maybe output list?
            }

            #region Implementation

            // todo - make private, send freqs over
            SymbolStats symbolStats = new();


            // Encoder symbol description
            // This (admittedly odd) selection of parameters was chosen to make
            // RansEncPutSymbol as cheap as possible.
            class CompressorSymbol
            {
                public uint x_max; // (Exclusive) upper bound of pre-normalization interval
                public uint rcp_freq; // Fixed-point reciprocal frequency
                public ushort rcp_shift; // Reciprocal shift
                public uint bias; // Bias
                public ushort cmpl_freq; // Complement of frequency: (1 << scale_bits) - freq
            }

            readonly CompressorSymbol[] symbols = new CompressorSymbol[256];


            // Initializes an encoder symbol to start "start" and frequency "freq"
            static void InitializeSymbol(CompressorSymbol symbol, uint cumulativeFrequency, uint frequency)
            {
                // sum of frequencies must be this
                const uint M = 1u << (int)ProbabilityBits;

                Debug.Assert(ProbabilityBits <= 16);
                Debug.Assert(cumulativeFrequency <= M);
                Debug.Assert(frequency <= M - cumulativeFrequency);

                // Say M := 1 << scale_bits.
                //
                // The original encoder does:
                //   x_new = (x/freq)*M + start + (x%freq)
                //
                // The fast encoder does (schematically):
                //   q     = mul_hi(x, rcp_freq) >> rcp_shift   (division)
                //   r     = x - q*freq                         (remainder)
                //   x_new = q*M + bias + r                     (new x)
                // plugging in r into x_new yields:
                //   x_new = bias + x + q*(M - freq)
                //        =: bias + x + q*cmpl_freq             (*)
                //
                // and we can just precompute cmpl_freq. Now we just need to
                // set up our parameters such that the original encoder and
                // the fast encoder agree.

                symbol.x_max = ((LowerBound >> (int)ProbabilityBits) << 8) * frequency;
                symbol.cmpl_freq = (ushort)(M - frequency);
                if (frequency < 2)
                {
                    // freq=0 symbols are never valid to encode, so it doesn't matter what
                    // we set our values to.
                    //
                    // freq=1 is tricky, since the reciprocal of 1 is 1; unfortunately,
                    // our fixed-point reciprocal approximation can only multiply by values
                    // smaller than 1.
                    //
                    // So we use the "next best thing": rcp_freq=0xffffffff, rcp_shift=0.
                    // This gives:
                    //   q = mul_hi(x, rcp_freq) >> rcp_shift
                    //     = mul_hi(x, (1<<32) - 1)) >> 0
                    //     = floor(x - x/(2^32))
                    //     = x - 1 if 1 <= x < 2^32
                    // and we know that x>0 (x=0 is never in a valid normalization interval).
                    //
                    // So we now need to choose the other parameters such that
                    //   x_new = x*M + start
                    // plug it in:
                    //     x*M + start                   (desired result)
                    //   = bias + x + q*cmpl_freq        (*)
                    //   = bias + x + (x - 1)*(M - 1)    (plug in q=x-1, cmpl_freq)
                    //   = bias + 1 + (x - 1)*M
                    //   = x*M + (bias + 1 - M)
                    //
                    // so we have start = bias + 1 - M, or equivalently
                    //   bias = start + M - 1.
                    symbol.rcp_freq = ~0u;
                    symbol.rcp_shift = 0;
                    symbol.bias = cumulativeFrequency + M - 1;
                }
                else
                {
                    // todo - test other format
                    // Alverson, "Integer Division using reciprocals"
                    // shift=ceil(log2(freq))
                    uint shift = 0;
                    while (frequency > 1u << (int)shift)
                        shift++;

                    symbol.rcp_freq = (uint)(((1UL << ((int)shift + 31)) + frequency - 1) / frequency);
                    symbol.rcp_shift = (ushort)(shift - 1);

                    // With these values, 'q' is the correct quotient, so we
                    // have bias=start.
                    symbol.bias = cumulativeFrequency;
                }
            }


            // Encodes the given symbol. 
            // NOTE: must encode in reverse order
            static void CompressSymbol(ref uint state, List<byte> outputBuffer, CompressorSymbol symbol)
            {

                // naively this function has a divide, but can remove by proper selection of reciprocals
                // renormalize
                // uint x = Renormalize(r, ptr, freq, scale_bits);
                // x = C(s,x)
                // r = ((x / freq) << (int)scale_bits) + (x % freq) + cumulativeFrequency;


                Debug.Assert(symbol.x_max != 0); // can't encode symbol with freq=0

                // renormalize the encoder
                while (state >= symbol.x_max)
                {
                    outputBuffer.Add((byte)state);
                    state >>= 8;
                }

                // x = C(s,x)
                // NOTE: written this way so we get a 32-bit "multiply high" when
                // available. If you're on a 64-bit platform with cheap multiplies
                // (e.g. x64), just bake the +32 into rcp_shift.
                uint q = (uint)(((ulong)state * symbol.rcp_freq) >> 32) >> symbol.rcp_shift;
                state = state + symbol.bias + q * symbol.cmpl_freq; // new state
            }

            static void Write4(uint value, ICollection<byte> compressedData)
            {
                compressedData.Add((byte)(value >> 24));
                compressedData.Add((byte)(value >> 16));
                compressedData.Add((byte)(value >> 8));
                compressedData.Add((byte)(value >> 0));

            }
            static void FlushCompressor(uint state, ICollection<byte> compressedData) => Write4(state, compressedData);

            #endregion
        }

        class Decompressor
        {

            public Decompressor()
            {
                for (var i = 0; i < 256; i++)
                    symbols[i] = new();
            }

            static uint Read4(IReadOnlyList<byte> data, ref int index)
            {
                uint val = data[index++];
                val |= (uint)(data[index++]) << 8;
                val |= (uint)(data[index++]) << 16;
                val |= (uint)(data[index++]) << 24;
                return val;
            }
            public byte[] Decompress(IReadOnlyList<byte> compressedData)
            {
                var index = 0; // data index

                var version = compressedData[index++]; // should be 0
                // need original size to know when to stop
                var decompressedLength = Read4(compressedData, ref index);
                var freqSize = compressedData[index++]; // should be 4

                uint[] scaledFrequencies = new uint[256];
                uint cum1 = 0, cum2 = 0; // track 2 of the cumulative frequency values for initialization

                // lookup cumulative value into symbol index
                // grows large as ProbabilityScaling grows
                // todo - other structures smaller and fast enough?
                byte[] cum2sym = new byte[ProbabilityScaling];

                for (var s = 0; s < 256; s++)
                {
                    // unpack 256 symbol frequencies for this block of data
                    var f = Read4(compressedData, ref index);
                    cum2 += f; // leading cumulative freq
                    scaledFrequencies[s] = f;

                    // cumulative->symbol table - todo - find faster?
                    for (var i = cum1; i < cum2; i++)
                        cum2sym[i] = (byte)s;

                    // NOTE: this uses the static model from the input data, needs shipped in compressed stream
                    InitializeSymbol(symbols[s], cum1, scaledFrequencies[s]);
                    
                    cum1 = cum2; // trailing cumulative freq
                }

                var decompressedBytes = new List<byte>();
                InitializeDecompressor(out var state, ref index, compressedData);

                for (var i = 0; i < decompressedLength; i++)
                {
                    var symbol = cum2sym[DecompressIndex(state)];
                    decompressedBytes.Add(symbol);
                    AdvanceState(ref state, compressedData, ref index, symbols[symbol]);
                }

                return decompressedBytes.ToArray(); // todo - slow, change up
            }

            // Decoder symbols are straightforward.
            class DecompressorSymbol
            {
                public ushort start; // Start of range.
                public ushort freq; // Symbol frequency.
            }
            DecompressorSymbol[] symbols = new DecompressorSymbol[256];

            // Initialize a decoder symbol to start "start" and frequency "freq"
            static void InitializeSymbol(DecompressorSymbol s, uint start, uint freq)
            {
                Debug.Assert(start <= (1 << 16));
                Debug.Assert(freq <= (1 << 16) - start);
                s.start = (ushort)start;
                s.freq = (ushort)freq;
            }


            // decoder
            // Initializes a rANS decoder.
            // Unlike the encoder, the decoder works forwards as you'd expect.
            static void InitializeDecompressor(out uint state, ref int index, IReadOnlyList<byte> data)
            {
                state = (uint)(data[index++] << 0);
                state |= (uint)(data[index++] << 8);
                state |= (uint)(data[index++] << 16);
                state |= (uint)(data[index++] << 24);
            }

            // Returns the current cumulative frequency (map it to a symbol yourself!)
            static uint DecompressIndex(uint r)
            {
                return r & ((1u << (int)ProbabilityBits) - 1);
            }

            // Advances in the bit stream by "popping" a single symbol with range start
            // "start" and frequency "freq". All frequencies are assumed to sum to "1 << scale_bits",
            // and the resulting bytes get written to ptr (which is updated).
            static void AdvanceState(ref uint state, IReadOnlyList<byte> data, ref int index, DecompressorSymbol sym)
            {
                uint start = sym.start;
                uint freq = sym.freq;
                uint mask = (1u << (int)ProbabilityBits) - 1;

                // s, x = D(x)
                state = freq * (state >> (int)ProbabilityBits) + (state & mask) - start;

                // renormalize
                while (state < LowerBound)
                {
                    state = (state << 8) | data[index];
                    ++index;
                }
            }
        }

        public static IList<byte> CompressBytes(IReadOnlyList<byte> data) => new Compressor().Compress(data);

        public static IList<byte> DecompressBytes(IReadOnlyList<byte> data) => new Decompressor().Decompress(data);
    }
}
