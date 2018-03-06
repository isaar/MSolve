using ISAAR.MSolve.Optimization.Algorithms.Metaheuristics.GeneticAlgorithms.Encodings;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Optimization.Commons
{
    public static class BinaryUtilities
    {
        /// <summary>
        /// Converts a binary subarray to a decimal integer 
        /// </summary>
        /// <param name="bits"></param>
        /// <param name="start"></param>
        /// <param name="end">Exclusive</param>
        /// <returns></returns>
        public static int StandardBinaryToDecimal(bool[] bits, int start, int length)
        {
            int dec = 0;
            for (int i = start; i < start + length; ++i)
            {
                dec += bits[i] ? dec + 1: dec;
                //Console.WriteLine(dec);
            }
            return dec;
        }

        public static void DecimalToStandardBinary(int dec, bool[] bits, int start, int length)
        {
            for (int i = start + length - 1; i >= start; --i)
            {
                bits[i] = (dec % 2 == 0) ? false : true;
                dec /= 2;
            }
        }

        /// <summary>
        /// Converts a Gray-coded binary subarray to a decimal integer 
        /// </summary>
        /// <param name="genotype"></param>
        /// <param name="start"></param>
        /// <param name="end">Exclusive</param>
        /// <returns></returns>
        public static int GrayCodeToDecimal(bool[] bits, int start, int length)
        {
            bool bin = bits[start];
            int dec = bin ? 1 : 0;
            for (int i = start + 1; i < start + length; ++i) // n-1 repetitions
            {
                bin = bin != bits[i];
                dec += bin ? dec + 1 : dec;
            }
            return dec;
        }

        public static void DecimalToGrayCode(int dec, bool[] bits, int start, int length)
        {
            bool[] binary = new bool[length];
            DecimalToStandardBinary(dec, binary, 0, length);
            // Binary to Gray code
            bits[start] = binary[0];
            for (int i = 1; i < length; ++i)
            {
                bits[start + i] = binary[i] != binary[i-1];
            }
        }

        public static void TestBinToDec()
        {
            bool[] bits1 = new bool[] { false, false, false };
            WriteDecimalRepresentations(bits1);
            bool[] bits2 = new bool[] { false, false, true };
            WriteDecimalRepresentations(bits2);
            bool[] bits3 = new bool[] { false, true, false };
            WriteDecimalRepresentations(bits3);
            bool[] bits4 = new bool[] { false, true, true };
            WriteDecimalRepresentations(bits4);
            bool[] bits5 = new bool[] { true, false, false };
            WriteDecimalRepresentations(bits5);
            bool[] bits6 = new bool[] { true, false, true };
            WriteDecimalRepresentations(bits6);
            bool[] bits7 = new bool[] { true, true, false };
            WriteDecimalRepresentations(bits7);
            bool[] bits8 = new bool[] { true, true, true };
            WriteDecimalRepresentations(bits8);
        }

        private static void WriteDecimalRepresentations(bool[] bits)
        {
            Console.Write("Bits = ");
            foreach (var entry in bits) Console.Write(entry ? 1 : 0);
            Console.Write(" -> from binary: " + StandardBinaryToDecimal(bits, 0, bits.Length));
            Console.WriteLine(" , from Gray code: " + GrayCodeToDecimal(bits, 0, bits.Length));
        }

        public static void TestDecToBin()
        {
            WriteBinaryRepresentations(0,3);
            WriteBinaryRepresentations(1,3);
            WriteBinaryRepresentations(2,3);
            WriteBinaryRepresentations(3,3);
            WriteBinaryRepresentations(4,3);
            WriteBinaryRepresentations(5,3);
            WriteBinaryRepresentations(6,3);
            WriteBinaryRepresentations(7,3);
            Console.WriteLine();
            Quantization quantization = new Quantization(3);
            quantization.PrintRanges();
            Console.WriteLine("\n double to int:");
            Console.WriteLine("0.55 -> " + quantization.NormalizedDoubleToInteger(0.55));
            Console.WriteLine("0.11 -> " + quantization.NormalizedDoubleToInteger(0.11));
            Console.WriteLine("0.95 -> " + quantization.NormalizedDoubleToInteger(0.95));
            Console.WriteLine("0.63 -> " + quantization.NormalizedDoubleToInteger(0.63));
        }

        private static void WriteBinaryRepresentations(int dec, int length)
        {
            bool[] binary = new bool[length];
            DecimalToStandardBinary(dec, binary, 0, length);
            bool[] grayCode = new bool[length];

            DecimalToGrayCode(dec, grayCode, 0, length);
            Console.Write("Decimal = " + dec);
            Console.Write(" -> Binary = ");
            foreach (var entry in binary) Console.Write(entry ? 1 : 0);
            Console.Write(" , Gray code = ");
            foreach (var entry in grayCode) Console.Write(entry ? 1 : 0);
            Console.WriteLine();
        }


    }
}
