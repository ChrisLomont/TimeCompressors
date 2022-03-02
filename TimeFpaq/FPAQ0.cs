// based on
// fpaq0 - Stationary order 0 file compressor.
// (C) 2004, Matt Mahoney under GPL, http://www.gnu.org/licenses/gpl.txt
// To compile: g++ -O fpaq0.cpp

// ported to C# to test against other compressors
// Chris Lomont 2022

using System.Diagnostics;

namespace TimeFpaq
{
    internal class FPAQ0
    {

        //////////////////////////// Predictor /////////////////////////

        /* A Predictor estimates the probability that the next bit of
           uncompressed data is 1.  Methods:
           p() returns P(1) as a 12 bit number (0-4095).
           update(y) trains the predictor with the actual bit (0 or 1).
        */

        class Predictor
        {
            int cxt = 0; // Context: last 0-8 bits with a leading 1
            int[,] ct = new int[512, 2]; // 0 and 1 counts in context cxt

            public Predictor()
            {
                cxt = 1;
            }

            // Assume a stationary order 0 stream of 9-bit symbols
            public int p()
            {
                return 4096 * (ct[cxt, 1] + 1) / (ct[cxt, 0] + ct[cxt, 1] + 2);
            }

            public void update(int y)
            {
                if (++ct[cxt, y] > 65534)
                {
                    ct[cxt, 0] >>= 1;
                    ct[cxt, 1] >>= 1;
                }

                if ((cxt += cxt + y) >= 512)
                    cxt = 1;
            }
        }


        /* An Encoder does arithmetic encoding.  Methods:
           Encoder(COMPRESS, f) creates encoder for compression to archive f, which
             must be open past any header for writing in binary mode
           Encoder(DECOMPRESS, f) creates encoder for decompression from archive f,
             which must be open past any header for reading in binary mode
           encode(bit) in COMPRESS mode compresses bit to file f.
           decode() in DECOMPRESS mode returns the next decompressed bit from file f.
           flush() should be called when there is no more to compress.
        */

        public enum Mode1
        {
            COMPRESS,
            DECOMPRESS
        }

        Mode1 Mode;


        Predictor predictor;
        Mode1 mode; // Compress or decompress?
        uint x1, x2; // Range, initially [0, 1), scaled by 2^32
        uint x; // Last 4 input bytes of archive.


        const int EOF = -1;

        // get a byte
        int getc()
        {
            if (pos >= input.Count) return EOF;
            return input[pos++];
        }

        // write byte out
        void putc(byte datum)
        {
            output.Add(datum);
            pos = output.Count;
        }


        // space, length of buffer
        public List<byte> output = new();
            public IReadOnlyList<byte> input;
        int pos = 0;

        // Constructor
        public FPAQ0(Mode1 m, IReadOnlyList<byte> input)
        {
            this.input = input;

            predictor = new();
            mode = m;
            x1 = 0;
            x2 = 0xffffffff;
            x = 0;

            // In DECOMPRESS mode, initialize x to the first 4 bytes of the archive
            if (mode == Mode1.DECOMPRESS)
            {
                for (int i = 0; i < 4; ++i)
                {
                    int c = getc();
                    if (c == EOF) c = 0;
                    x = (x << 8) + (uint)(c & 0xff);
                }
            }
        }

        public void encodeByte(byte b)
        {

            encode(0); // mark as data byte
            for (int i = 7; i >= 0; --i)
                encode((b >> i) & 1);

        }

        /* encode(y) -- Encode bit y by splitting the range [x1, x2] in proportion
        to P(1) and P(0) as given by the predictor and narrowing to the appropriate
        subrange.  Output leading bytes of the range as they become known. */

        public void encode(int y)
        {

            // Update the range
            uint xmid = (uint)(x1 + ((x2 - x1) >> 12) * predictor.p());
            Trace.Assert(xmid >= x1 && xmid < x2);
            if (y != 0)
                x2 = xmid;
            else
                x1 = xmid + 1;
            predictor.update(y);

            // Shift equal MSB's out
            while (((x1 ^ x2) & 0xff000000) == 0)
            {
                putc((byte)(x2 >> 24));
                x1 <<= 8;
                x2 = (x2 << 8) + 255;
            }
        }

        bool done = false;
        // decode byte, or -1
        public int DecodeByte()
        {
            var t = decode();
            if (t == 0 || done)
            {
                done = true;
                return -1;
            }

            // get a byte
            int c = 1;
            while (c < 256)
                c += c + decode();
            return c - 256;
        }

        /* Decode one bit from the archive, splitting [x1, x2] as in the encoder
        and returning 1 or 0 depending on which subrange the archive point x is in.
        */
        public int decode()
        {

            // Update the range
            uint xmid = (uint)(x1 + ((x2 - x1) >> 12) * predictor.p());
            Trace.Assert(xmid >= x1 && xmid < x2);
            int y = 0;
            if (x <= xmid)
            {
                y = 1;
                x2 = xmid;
            }
            else
                x1 = xmid + 1;

            predictor.update(y);

            // Shift equal MSB's out
            while (((x1 ^ x2) & 0xff000000) == 0)
            {
                x1 <<= 8;
                x2 = (x2 << 8) + 255;
                int c = getc();
                if (c == EOF) c = 0;
                x = (x << 8) + (uint)c;
            }

            return y;
        }

        // Should be called when there is no more to compress
        public void flush()
        {

            // In COMPRESS mode, write out the remaining bytes of x, x1 < x < x2
            if (mode == Mode1.COMPRESS)
            {
                encode(1); // mark end of data bytes

                while (((x1 ^ x2) & 0xff000000) == 0)
                {
                    putc((byte)(x2 >> 24));
                    x1 <<= 8;
                    x2 = (x2 << 8) + 255;
                }

                putc((byte)(x2 >> 24)); // First unequal byte
            }
        }
    }

}

