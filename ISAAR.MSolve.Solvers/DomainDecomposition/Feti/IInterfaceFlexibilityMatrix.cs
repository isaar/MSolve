﻿using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.LinearAlgebra.Vectors;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Feti
{
    internal interface IInterfaceFlexibilityMatrix
    {
        int Order { get; }
        void Multiply(Vector lhs, Vector rhs);
    }
}
