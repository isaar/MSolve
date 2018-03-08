using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.Materials.Interfaces
{
    public interface IFiniteElementMaterial : ICloneable
    {
        int ID { get; }
        bool Modified { get; }
        void ResetModified();
        double[] Coordinates { get; set; }
    }
}
