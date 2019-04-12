using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public enum DOFType
    {
        Unknown = 0,
        X = 1,
        Y = 2,
        Z = 3,
        RotX = 4,
        RotY = 5,
        RotZ = 6,
        Pore = 7,
        Temperature = 8
    }
}
