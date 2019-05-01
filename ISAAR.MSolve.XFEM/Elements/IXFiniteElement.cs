using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Elements
{
    public interface IXFiniteElement : IElement, IElementType
    {
        IReadOnlyList<XNode> Nodes { get; }
        XSubdomain Subdomain { get; set; }
    }
}