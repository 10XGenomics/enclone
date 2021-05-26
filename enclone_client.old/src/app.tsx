import Ansi from 'ansi-to-react';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faSpinner } from '@fortawesome/free-solid-svg-icons';
import * as React from 'react';
import * as ReactDOM from 'react-dom';

import Plot from './components/Plot';

import { EncloneApi } from './preload'
declare global {
  interface Window { encloneApi: EncloneApi; }
}

export const App = (): JSX.Element => {
  const [state, setState] = React.useState({
    argsLoaded: '', // loaded enclone args
    clonotypeId: '', // searched for clonotype id
    encloneArgs: '', // user-inputted enclone args
    encloneError: false,
    encloneLoading: false,
    plot: '', // returned svg data
    table: '',
  });

  function handleEncloneReq() {
    if (!state.encloneArgs) {
      return;
    }
    setState({
      ...state,
      argsLoaded: '',
      clonotypeId: '',
      encloneError: false,
      encloneLoading: true,
      plot: '',
      table: '',
    });
    window.encloneApi.enclone(state.encloneArgs);
  }

  function handleGetClonotype() {
    if (state.argsLoaded) {
      window.encloneApi.getClonotypesById(state.clonotypeId);
    }
  }

  // Create the process event listener once on load
  React.useEffect(() => {
    window.encloneApi.onClonotypesByIdResponse((data: any) => {
      if (!data || !data.table) {
        setState(state => ({ ...state, table: '' }));
        return;
      }
      setState(state => ({ ...state, table: data.table }));
    });
    window.encloneApi.onEncloneResponse((data: any) => {
      if (data instanceof Error) {
        setState(state => ({
          ...state,
          argsLoaded: '',
          encloneArgs: '',
          encloneError: true,
          encloneLoading: false,
        }));
        return;
      }
      setState(state => ({
        ...state,
        argsLoaded: data.args,
        encloneError: false,
        encloneLoading: false,
        plot: data.plot,
        table: data.table,
      }));
    });
  }, []);

  return (
    <div>
      <h3>enclone GUI</h3>
      <div className="instruction">
        <span>Run enclone with args first:</span>
        <input
          onChange={e => setState({ ...state, encloneArgs: e.target.value })}
          placeholder="Type enclone params"
          style={{ width: '260px' }}
          type="text"
          value={state.encloneArgs}
        />
        <button disabled={state.encloneLoading} onClick={handleEncloneReq}>
          Send
        </button>
        {state.encloneLoading && <FontAwesomeIcon className="fa-spin" icon={faSpinner} />}
        <div className="hint">
          If you're not sure what to type, try: BCR=~/enclone/datasets/1031851
        </div>
      </div>
      {state.encloneError && (
        <div className="warning">
          Something went wrong. If the error persists, restart the server
        </div>
      )}
      {state.argsLoaded && (
        <>
          <div className="success">{`Enclone running with args ${state.argsLoaded}`}</div>
          <div className="instruction" style={{ overflowY: 'scroll' }}>
            <span>Clonotype Idx: </span>
            <input
              onChange={e => setState({ ...state, clonotypeId: e.target.value })}
              placeholder="Type a Clonotype Id"
              type="text"
              value={state.clonotypeId}
            />
            <button onClick={handleGetClonotype}>Send</button>
          </div>
          <div className="data-container">
            <Plot svg={state.plot} />
          </div>
          <div className="data-container">
            <Ansi>{state.table}</Ansi>
          </div>
        </>
      )}
    </div>
  );
};

const rootElement = document.getElementById('root');
ReactDOM.render(<App />, rootElement);
