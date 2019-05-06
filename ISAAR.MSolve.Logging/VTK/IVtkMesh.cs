using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Logging.VTK
{
    public interface IVtkMesh
    {
        IReadOnlyList<VtkCell> VtkCells { get; }
        IReadOnlyList<VtkPoint> VtkPoints { get; }
    }
}
