// test FPAQ0 versus my arithmetic compressor
// looks like mine 25% faster for files...., slightly worse compression (fix that!)
// rANS about 10 times faster! figure out how to use in general :)

using System.Diagnostics;
using TimeFpaq;

//TestDiv.RunTests();
//return;

var compressors = new Type[]
{
    typeof(Fpaq0), // decompressor broken :(
    typeof(Arith), 
    typeof(rAns2)
};

Run(testDecompressor:true);
return;

void Run(bool testDecompressor = false)
{
    var states = new List<State>();
    foreach (var c in compressors)
        states.Add(new State { type = c });

    long totalBytes = 0, errors = 0;
foreach (var f in Directory.EnumerateFiles(@"e:\662BinariesOLD", "*.dll"))
//foreach (var f in Directory.EnumerateFiles(@"i:\671\software","*.*", SearchOption.AllDirectories))
//    foreach (var f in Directory.EnumerateFiles(@"..\..\..\..\", "*.*", SearchOption.AllDirectories))
    //foreach (var f in Directory.EnumerateFiles(@"d:\", "*.*", SearchOption.AllDirectories))
    {
        byte[] data = new byte[0];
        bool error = false;
        try
        {
            if (new FileInfo(f).Length > 100_000_000)
            {
                continue;
            }

            data = File.ReadAllBytes(f);
            Console.WriteLine($"File {Path.GetFileName(f)}: {data.Length} bytes: ({errors} errors) ");

            foreach (var s in states)
            {
                var inst = (ICompressor)Activator.CreateInstance(s.type);
                var sw = new Stopwatch();
                sw.Start();
                var compressedData = inst.Compress(data);
                sw.Stop();
                s.compressionElapsedTimeSingle = sw.Elapsed;
                s.compressedBytesSingle = compressedData.Count;
                if (testDecompressor)
                {
                    // in case does not play well
                    inst = (ICompressor)Activator.CreateInstance(s.type);

                    IReadOnlyList<byte> rd = ((List<byte>)compressedData).AsReadOnly();
                    sw.Restart();
                    var decompressedData = inst.Decompress(rd);
                    sw.Stop();
                    s.decompressionElapsedTimeSingle = sw.Elapsed;
                    s.decompressedBytesSingle = decompressedData.Count;

                    if (!Matches(decompressedData, data))
                    {
                        throw new Exception($"mismatch error in {s.type.Name}!");
                    }
                }
            }
        }
        catch (Exception ex)
        {
            Console.WriteLine($"Exception {ex}");
            error = true;
            errors++;
        }

        if (!error)
        {
            totalBytes += data.Length;

            foreach (var s in states)
            {
                s.compressedBytesTotal += s.compressedBytesSingle;
                s.decompressedBytesTotal += s.decompressedBytesSingle;
                s.compressionElapsedTime += s.compressionElapsedTimeSingle;
                s.decompressionElapsedTime += s.decompressionElapsedTimeSingle;
                Console.Write(
                    $"   {s.type.Name} : {s.compressedBytesSingle}, {totalBytes * 1000 / s.compressionElapsedTime.Milliseconds:D9} bytes/sec, {s.compressedBytesTotal * 100 / totalBytes}% ratio");
                if (testDecompressor)
                {
                    Console.Write(
                        $" , dec {totalBytes * 1000 / s.decompressionElapsedTime.Milliseconds:D9} bytes/sec");
                }

                Console.WriteLine();

            }
        }


        Console.WriteLine();
    }
}

bool Matches(IList<byte> data1, IList<byte> data2)
{
    if (data1.Count != data2.Count)
        return false;
    for (var i = 0; i < data1.Count; ++i)
        if (data1[i] != data2[i])
            return false;

    return true;
}

class State
{
    public TimeSpan compressionElapsedTime = TimeSpan.Zero, decompressionElapsedTime = TimeSpan.Zero;
    public TimeSpan compressionElapsedTimeSingle = TimeSpan.Zero, decompressionElapsedTimeSingle = TimeSpan.Zero;
    public long compressedBytesTotal, decompressedBytesTotal;
    public long compressedBytesSingle, decompressedBytesSingle;
    public Type type;
}

class Arith : ICompressor
{
    public IList<byte> Compress(IReadOnlyList<byte> data) => Lomont.Compression.Compressor.CompressBytes(data);

    public IList<byte> Decompress(IReadOnlyList<byte> data) => Lomont.Compression.Decompressor.DecompressBytes(data);
}

class Fpaq0 : ICompressor
{
    public IList<byte> Compress(IReadOnlyList<byte> data)
    {
        var fp = new FPAQ0(FPAQ0.Mode1.COMPRESS, data);
        foreach (var b in data)
            fp.encodeByte(b);
        fp.flush();
        return fp.output;
    }

    public IList<byte> Decompress(IReadOnlyList<byte> data)
    {

        var de = new FPAQ0(FPAQ0.Mode1.DECOMPRESS, data);

        var output = new List<byte>();
        while (true)
        {
            var c = de.DecodeByte();
            if (c == -1)
                break;
            output.Add((byte)c);
        }

        return output;
    }
}

class rAns2 : ICompressor
{
    public IList<byte> Compress(IReadOnlyList<byte> data) => rANS.CompressBytes(data);

    public IList<byte> Decompress(IReadOnlyList<byte> data) => rANS.DecompressBytes(data);
}

interface ICompressor
{

    // compress data
    IList<byte> Compress(IReadOnlyList<byte> data);
    IList<byte> Decompress(IReadOnlyList<byte> data);
}

