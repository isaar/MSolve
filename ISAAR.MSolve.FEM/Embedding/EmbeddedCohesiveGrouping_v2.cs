using ISAAR.MSolve.FEM.Embedding;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;

//TODO: this and EmbeddedGrouping have most things in common. Use a base class for them and template method or use polymorhism from the composed classes.
namespace ISAAR.MSolve.FEM.Embedding
{
    public class EmbeddedCohesiveGrouping_v2
    {
        private readonly Model_v2 model;
        private readonly IEnumerable<Element> hostGroup;
        private readonly IEnumerable<Element> embeddedGroup;
        private readonly bool hasEmbeddedRotations = false;

        public IEnumerable<Element> HostGroup { get { return hostGroup; } }
        public IEnumerable<Element> EmbeddedGroup { get { return embeddedGroup; } }

        public EmbeddedCohesiveGrouping_v2(Model_v2 model, IEnumerable<Element> hostGroup, IEnumerable<Element> embeddedGroup, bool hasEmbeddedRotations)
        {
            this.model = model;
            this.hostGroup = hostGroup;
            this.embeddedGroup = embeddedGroup;
            this.hasEmbeddedRotations = hasEmbeddedRotations;
            hostGroup.Select(e => e.ElementType).Distinct().ToList().ForEach(et =>
            {
                if (!(et is IEmbeddedHostElement))
                    throw new ArgumentException("EmbeddedGrouping: One or more elements of host group does NOT implement IEmbeddedHostElement.");
            });
            embeddedGroup.Select(e => e.ElementType).Distinct().ToList().ForEach(et =>
            {
                if (!(et is IEmbeddedElement))
                    throw new ArgumentException("EmbeddedGrouping: One or more elements of embedded group does NOT implement IEmbeddedElement.");
            });
            UpdateNodesBelongingToEmbeddedElements();
        }

        public EmbeddedCohesiveGrouping_v2(Model_v2 model, IEnumerable<Element> hostGroup, IEnumerable<Element> embeddedGroup)
            : this(model, hostGroup, embeddedGroup, false)
        {
        }

        private void UpdateNodesBelongingToEmbeddedElements()
        {
            IEmbeddedDOFInHostTransformationVector transformer;
            if (hasEmbeddedRotations)
                //transformer = new Hexa8TranslationAndRotationTransformationVector();
                throw new NotImplementedException();
            else
                transformer = new Hexa8LAndNLTranslationTransformationVector();

            foreach (var embeddedElement in embeddedGroup)
            {
                var elType = (IEmbeddedElement)embeddedElement.ElementType;
                foreach (var node in embeddedElement.Nodes.Skip(8))
                {
                    var embeddedNodes = hostGroup
                        .Select(e => ((IEmbeddedHostElement)e.ElementType).BuildHostElementEmbeddedNode(e, node, transformer))
                        .Where(e => e != null);
                    foreach (var embeddedNode in embeddedNodes)
                    {
                        if (elType.EmbeddedNodes.Count(x => x.Node == embeddedNode.Node) == 0)
                            elType.EmbeddedNodes.Add(embeddedNode);

                        // Update embedded node information for elements that are not inside the embedded group but contain an embedded node.
                        foreach (var element in model.Elements.Except(embeddedGroup))
                            if (element.ElementType is IEmbeddedElement && element.Nodes.Contains(embeddedNode.Node))
                            {
                                var currentElementType = (IEmbeddedElement)element.ElementType;
                                if (!currentElementType.EmbeddedNodes.Contains(embeddedNode))
                                {
                                    currentElementType.EmbeddedNodes.Add(embeddedNode);
                                    element.ElementType.DOFEnumerator = new CohesiveElementEmbedder_v2(model, element, transformer);
                                }
                            }
                    }
                }

                embeddedElement.ElementType.DOFEnumerator = new CohesiveElementEmbedder_v2(model, embeddedElement, transformer);
            }
        }
    }
}
