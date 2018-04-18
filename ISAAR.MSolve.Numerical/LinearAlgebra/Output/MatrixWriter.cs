using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Output
{
    // TODO: properties regarding precision, justification (perhaps in concrete) etc.
    public abstract class MatrixWriter
    {
        public void WriteToConsole()
        {
            // Point the stream to standard output
            var writer = new StreamWriter(Console.OpenStandardOutput());
            writer.AutoFlush = true;
            Console.SetOut(writer);

            // Call the abstract method
            WriteToStream(writer);

            // Recover the standard output stream
            var standardOutput = new StreamWriter(Console.OpenStandardOutput());
            standardOutput.AutoFlush = true;
            Console.SetOut(standardOutput);
        }

        /// <summary>
        /// Write the entries of the matrix/vector/array to a specified file. If the file doesn't exist a new one will be created.
        /// </summary>
        /// <param name="path">The path of the file and its extension.</param>
        /// <param name="append">If the file already exists: Pass <see cref="append"/> = true to write after the current end of 
        ///     the file. Pass<see cref="append"/> = false to overwrite the file.</param>
        public void WriteToFile(string path, bool append = false)
        {
            using (var writer = new StreamWriter(path, append))
            {
#if DEBUG
                writer.AutoFlush = true; // To look at intermediate output at certain breakpoints
#endif
                WriteToStream(writer);
            }
        }

        

        protected abstract void WriteToStream(StreamWriter writer);
    }
}
