using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.Numerical.LinearAlgebra.Logging
{
    // TODO: Perhaps writing/reading matrices to/from files/console must be done by specialized classes, not inside the matrix 
    // classes.
    public interface IWriteable
    {
        void WriteToConsole();
        void WriteToFile(string path, bool append = false);
    }
}
