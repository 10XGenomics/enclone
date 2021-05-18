import * as React from 'react';

interface PlotProps {
  svg: string;
}

const Plot = ({ svg }: PlotProps): JSX.Element => {
  const [tooltipText, setTooltipText] = React.useState('{}');

  const svgWrapper = React.useRef(null);

  React.useEffect(() => {
    if (svgWrapper.current) {
      const svgNode = svgWrapper.current.getElementsByTagName('svg')[0];
      if (svgNode) {
        svgNode.addEventListener('mouseover', (event: MouseEvent) => {
          const target = event.target as SVGElement;
          const tooltip = target.dataset && target.dataset.tooltip;
          if (tooltip) {
            setTooltipText(tooltip);
          }
        });
      }
    }
  }, [svg]);

  const tooltipEntries = Object.entries(JSON.parse(tooltipText));

  if (svg && svg.length > 0) {
    return (
      <div className="viz">
        <div dangerouslySetInnerHTML={{ __html: svg }} ref={svgWrapper} />
        <div>
          {Boolean(tooltipEntries.length) && (
            <div className="tooltip-box">
              {tooltipEntries.map(([k, v]) => {
                return <div key={k}>{`${k}: ${v}`}</div>;
              })}
            </div>
          )}
        </div>
      </div>
    );
  }

  return null;
};

export default Plot;
