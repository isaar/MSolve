using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.Discretization.Interfaces
{
    public interface IElement
    {
		int ID { get; set; }
		ISubdomain Subdomain { get; set; }
		int[] DOFs { get; }
		IList<INode> Nodes { get; }
	    Dictionary<int, INode> NodesDictionary { get; }
		IElementType ElementType { get; set; }
	}
}
