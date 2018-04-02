using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Logging
{
    public interface IWriteable
    {
        void WriteToConsole();
        void WriteToFile(string path, bool append = false);
    }
}
