using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IFiniteElementMaterial
    {
        int ID { get; }
        bool Modified { get; }
        void ResetModified();
        double[] Coordinates { get; set; }
        double YoungModulus { get; }
        double PoissonRatio { get; }

        void SaveState();
        void ClearState();
        void ClearStresses();
    }
}
