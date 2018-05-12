using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IElement
    {
	    int ID { get; set; }
		IElementType IElementType { get; }
	    IList<INode> INodes { get; }
    }
}
