using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Embedding
{
    public class EmbeddedGrouping_v2
    {
        private readonly Model_v2 model;
        private readonly IEnumerable<Element_v2> hostGroup;
        private readonly IEnumerable<Element_v2> embeddedGroup;
        private readonly bool hasEmbeddedRotations = false;

        public IEnumerable<Element_v2> HostGroup { get { return hostGroup; } }
        public IEnumerable<Element_v2> EmbeddedGroup { get { return embeddedGroup; } }

        public EmbeddedGrouping_v2(Model_v2 model, IEnumerable<Element_v2> hostGroup, IEnumerable<Element_v2> embeddedGroup, bool hasEmbeddedRotations)
        {
            this.model = model;
            this.hostGroup = hostGroup;
            this.embeddedGroup = embeddedGroup;
            this.hasEmbeddedRotations = hasEmbeddedRotations;
            hostGroup.Select(e => e.ElementType).Distinct().ToList().ForEach(et =>
            {
                if (!(et is IEmbeddedHostElement_v2))
                    throw new ArgumentException("EmbeddedGrouping: One or more elements of host group does NOT implement IEmbeddedHostElement.");
            });
            embeddedGroup.Select(e => e.ElementType).Distinct().ToList().ForEach(et =>
            {
                if (!(et is IEmbeddedElement_v2))
                    throw new ArgumentException("EmbeddedGrouping: One or more elements of embedded group does NOT implement IEmbeddedElement.");
            });
            UpdateNodesBelongingToEmbeddedElements();
        }

        public EmbeddedGrouping_v2(Model_v2 model, IEnumerable<Element_v2> hostGroup, IEnumerable<Element_v2> embeddedGroup)
            : this(model, hostGroup, embeddedGroup, false)
        {
        }

        private void UpdateNodesBelongingToEmbeddedElements()
        {
            IEmbeddedDOFInHostTransformationVector_v2 transformer;
            if (hasEmbeddedRotations)
                transformer = new Hexa8TranslationAndRotationTransformationVector_v2();
            else
                transformer = new Hexa8LAndNLTranslationTransformationVector_v2();

            foreach (var embeddedElement in embeddedGroup)
            {
                var elType = (IEmbeddedElement_v2)embeddedElement.ElementType;
                foreach (var node in embeddedElement.Nodes)
                {
                    var embeddedNodes = hostGroup
                        .Select(e => ((IEmbeddedHostElement_v2)e.ElementType).BuildHostElementEmbeddedNode(e, node, transformer))
                        .Where(e => e != null);
                    foreach (var embeddedNode in embeddedNodes)
                    {
                        if (elType.EmbeddedNodes.Count(x => x.Node == embeddedNode.Node) == 0)
                            elType.EmbeddedNodes.Add(embeddedNode);

                        // Update embedded node information for elements that are not inside the embedded group but contain an embedded node.
                        foreach (var element in model.Elements.Except(embeddedGroup))
                            if (element.ElementType is IEmbeddedElement_v2 && element.Nodes.Contains(embeddedNode.Node))
                            {
                                var currentElementType = (IEmbeddedElement_v2)element.ElementType;
                                if (!currentElementType.EmbeddedNodes.Contains(embeddedNode))
                                {
                                    currentElementType.EmbeddedNodes.Add(embeddedNode);
                                    element.ElementType.DofEnumerator = new ElementEmbedder_v2(model, element, transformer);
                                }
                            }
                    }
                }

                embeddedElement.ElementType.DofEnumerator = new ElementEmbedder_v2(model, embeddedElement, transformer);
            }
        }
    }
}
