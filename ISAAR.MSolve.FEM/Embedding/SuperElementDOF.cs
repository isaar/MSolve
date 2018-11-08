using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using System;

namespace ISAAR.MSolve.FEM.Embedding
{
    public class SuperElementDOF
    {
        public Element Element { get; set; }
        public INode HostNode { get; set; }
        public INode EmbeddedNode { get; set; }
        public DOFType DOF { get; set; }

        public override bool Equals(object obj)
        {
            if (obj is SuperElementDOF == false) return false;
            var e = obj as SuperElementDOF;
            if (e == null) return false;

            if ((e.Element == null && this.Element != null) || (e.Element != null && this.Element == null)) return false;
            if ((e.HostNode == null && this.HostNode != null) || (e.HostNode != null && this.HostNode == null)) return false;
            if ((e.EmbeddedNode == null && this.EmbeddedNode != null) || (e.EmbeddedNode != null && this.EmbeddedNode == null)) return false;

            return (e.DOF == this.DOF && 
                ((e.Element == null && this.Element == null) || (e.Element.ID == this.Element.ID)) && 
                ((e.HostNode == null && this.HostNode == null) || (e.HostNode.ID == this.HostNode.ID)) &&
                ((e.EmbeddedNode == null && this.EmbeddedNode == null) || (e.EmbeddedNode.ID == this.EmbeddedNode.ID)));
        }

        public override int GetHashCode()
        {
            int elementID = Element == null ? 0 : Element.ID;
            int hostNodeID = HostNode == null ? 0 : HostNode.ID;
            return String.Format("H{3}E{0}N{1}{2}", elementID, hostNodeID, DOF, EmbeddedNode.ID).GetHashCode();
        }

        public override string ToString()
        {
            int elementID = Element == null ? 0 : Element.ID;
            int hostNodeID = HostNode == null ? 0 : HostNode.ID;
            return String.Format("N:{3} -> E:{0}, N:{1}, {2}", elementID, hostNodeID, DOF, EmbeddedNode.ID);
        }
    }

}
