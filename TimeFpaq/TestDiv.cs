using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TimeFpaq
{
    // test some reciprocal division routines
    static class TestDiv
    {
        public static void RunTests()
        {
            //var x1 = uint.MaxValue-2;
            //var d1 = uint.MaxValue-1; // fail case on Div1, works Div2
            //Test1(x1,d1,Div2); // seems to work (prove all via Z3?)
            //return;

            TestBatch(10000, Div1, false); // fails with test high true (off by 1 at most?)
            TestBatch(10000, Div2, true);

            void TestBatch(int passes, Func<uint, (uint, uint, int)> div, bool testHigh)
            {
                var r = new Random(1234);
                int maxB = 21; // bits - test 0 to 1<<w
                for (var pass = 0; pass < passes; ++pass)
                {
                    var d = (uint)(r.Next() + 1);

                    if (r.Next(10) > 5)
                        d = (uint)pass + 1; // some sequential

                    if (testHigh && r.Next(10) > 5)
                        d = uint.MaxValue - d;

                    var (a, b, s) = div(d); // or Div1
                    for (var i = 0UL; i < 1UL << maxB; ++i)
                    {
                        var x = (uint)i;
                        if (testHigh && r.Next(10) > 5)
                            x = uint.MaxValue - x;
                        if (!Test(x, d, a, b, s))
                            break;
                    }

                    if ((pass % (passes / 50)) == 0)
                        Console.WriteLine($"Pass {pass + 1}/{passes}");
                }
            }

            void Test1(uint x, uint d, Func<uint, (uint, uint, int)> div)
            {
                var (a, b, s) = div(d);
                Test(x, d, a, b, s);
            }

            bool Test(uint x, uint d, uint a, uint b, int s)
            {
                //x = (uint)r.Next();
                var (err, q1, q2) = Compute(x, d, a, b, s);
                if (err > 0)
                {
                    Console.WriteLine($"Div {d} => {a},{b},{s}");
                    Console.WriteLine($"Error: {x}/{d} => {q1} != {q2}");
                    return false;
                }

                return true;
            }

            (long err, uint q1, uint q2) Compute(uint x, uint d, uint a, uint b, int s)
            {
                //x = (uint)r.Next();
                var q1 = (uint)(x / d);
                var q2 = (uint)((a * (ulong)x + b) >> s);
                return (Math.Abs((long)q1 - (long)q2), q1, q2);
            }


            // Cavagnino - 2007 - Efficient Algorithms for Integer Division by Constants Using Multiplication

            // compute floor(x/d) for x in [0,2^N-1], fixed d, using form
            // Floor[(a*x+b)/2^s]
            // given d, compute a,b,s
            (uint a, uint b, int s) Div3(uint d)
            {
                // Moller 2010 - "Improved division by invariant integers"

                // L = word length
                // b=2^L word size
                //
                // quotient
                // q=Floor[U/d], r = U-qd
                // 
                // assumes d has most sig bit set - todo - shift in?
                // b/2 <= d < b  => 1/b < 1/d <= 2/b
                //
                // let v = floor((b*b-1)/d) - b
                // then

                var L = 32;
                UInt64 b = 1UL << L;
                var v = UInt64.MaxValue / d - b;

                // seems to be only for double word divided by word


                return (0, 0, 0);
            }

            // Theo Drane (2012) - compute using minimal values? no other benefit?


            // Granlund: gives result of form
            // d=> (a,s1,s2)
            // t1 = (a*x)>>32;
            // q2 =  (t1 + (x-t1)>>s1)>>s2
            // maybe removes the 2w+1 bit requrements?



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

            }



            // compute floor(x/d) for x in [0,2^N-1], fixed d, using form
            // Floor[(a*x+b)/2^s]
            // given d, compute a,b,s
            // NOTE: b = 0 always in this method, suffers from needing one more bit in reciprocal than in possible numerator or divisor
            (uint a, uint b, int s) Div1(uint y)
            {
                // Alverson, 1991 "Integer Division using reciprocals"
                // w = wordlength
                // s = w+ceil(log2(y))
                // a = ceil(2^s/y)
                // requires w+1 bit accuracy for a
                // so, to make fit in 32 bit words, restrict numerator to 31 bits, 
                // reduce shift by one

                int shift = 0;
                if (y > (1U << 31))
                    shift = 32;
                else
                {
                    while (y > 1u << shift)
                        shift++;
                }

                var a = (uint)(((1UL << ((int)shift + 31)) + y - 1) / y); // ceiling,compute as 64 bit division
                return (a, 0, shift - 1 + 32);
            }


        }
    }
}
