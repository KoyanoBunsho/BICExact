import { ChangeEvent, FormEvent, useMemo, useState } from 'react';
import {
  Activity,
  ArrowRight,
  CheckCircle2,
  FileUp,
  FlaskConical,
  Loader2,
  RotateCcw,
} from 'lucide-react';

type Criterion = 'bic' | 'aic';
type Method = 'exact' | 'lh';

type EstimateResult = {
  results: Record<string, string>;
  rawOutput?: string;
};

type FormState = {
  file1: File | null;
  file2: File | null;
  chain1: string;
  chain2: string;
  criterion: Criterion;
  method: Method;
};

const initialForm: FormState = {
  file1: null,
  file2: null,
  chain1: 'A',
  chain2: 'A',
  criterion: 'bic',
  method: 'exact',
};

const examples = [
  { label: 'Criterion', value: 'BIC / AIC' },
  { label: 'Method', value: 'Exact / LH' },
  { label: 'Input', value: 'Two PDB conformations' },
];

function App() {
  const [form, setForm] = useState<FormState>(initialForm);
  const [result, setResult] = useState<EstimateResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [isSubmitting, setIsSubmitting] = useState(false);

  const canSubmit = useMemo(
    () =>
      Boolean(form.file1 && form.file2 && form.chain1.trim() && form.chain2.trim()) &&
      !isSubmitting,
    [form, isSubmitting],
  );

  const setFile = (key: 'file1' | 'file2') => (event: ChangeEvent<HTMLInputElement>) => {
    setForm((current) => ({ ...current, [key]: event.target.files?.[0] ?? null }));
  };

  const submit = async (event: FormEvent<HTMLFormElement>) => {
    event.preventDefault();
    if (!form.file1 || !form.file2) return;

    setError(null);
    setResult(null);
    setIsSubmitting(true);

    const body = new FormData();
    body.append('file1', form.file1);
    body.append('file2', form.file2);
    body.append('chain1', form.chain1.trim());
    body.append('chain2', form.chain2.trim());
    body.append('bic', form.criterion);
    body.append('exactmo', form.method);

    try {
      const response = await fetch('/api/estimate', {
        method: 'POST',
        body,
      });
      const data = await response.json().catch(() => null);

      if (!response.ok) {
        throw new Error(data?.error ?? `Request failed with status ${response.status}`);
      }

      setResult(data as EstimateResult);
    } catch (caught) {
      setError(caught instanceof Error ? caught.message : 'Estimation failed.');
    } finally {
      setIsSubmitting(false);
    }
  };

  const reset = () => {
    setForm(initialForm);
    setResult(null);
    setError(null);
  };

  return (
    <main className="app-shell">
      <section className="hero">
        <div className="hero-copy">
          <span className="eyebrow">
            <Activity size={16} aria-hidden="true" />
            Protein hinge estimation
          </span>
          <h1>BIC Exact</h1>
          <p>
            Estimate hinge counts and positions from two conformations of the same protein.
          </p>
        </div>

        <div className="hero-metrics" aria-label="Application summary">
          {examples.map((item) => (
            <div className="metric" key={item.label}>
              <span>{item.label}</span>
              <strong>{item.value}</strong>
            </div>
          ))}
        </div>
      </section>

      <section className="workspace" aria-label="Estimator">
        <form className="panel estimator" onSubmit={submit}>
          <div className="panel-header">
            <div>
              <span className="section-kicker">Run estimation</span>
              <h2>Upload structures</h2>
            </div>
            <button className="icon-button" type="button" onClick={reset} title="Reset form">
              <RotateCcw size={18} aria-hidden="true" />
            </button>
          </div>

          <div className="upload-grid">
            <FileInput
              id="file1"
              label="PDB file 1"
              file={form.file1}
              onChange={setFile('file1')}
            />
            <FileInput
              id="file2"
              label="PDB file 2"
              file={form.file2}
              onChange={setFile('file2')}
            />
          </div>

          <div className="field-grid">
            <label>
              <span>Chain ID 1</span>
              <input
                value={form.chain1}
                maxLength={8}
                onChange={(event) => setForm({ ...form, chain1: event.target.value })}
                required
              />
            </label>
            <label>
              <span>Chain ID 2</span>
              <input
                value={form.chain2}
                maxLength={8}
                onChange={(event) => setForm({ ...form, chain2: event.target.value })}
                required
              />
            </label>
          </div>

          <div className="control-row" aria-label="Information criterion">
            <span>Criterion</span>
            <SegmentedControl
              value={form.criterion}
              options={[
                { label: 'BIC', value: 'bic' },
                { label: 'AIC', value: 'aic' },
              ]}
              onChange={(criterion) => setForm({ ...form, criterion })}
            />
          </div>

          <div className="control-row" aria-label="Estimation method">
            <span>Method</span>
            <SegmentedControl
              value={form.method}
              options={[
                { label: 'Exact', value: 'exact' },
                { label: 'LH', value: 'lh' },
              ]}
              onChange={(method) => setForm({ ...form, method })}
            />
          </div>

          {error && <div className="message error">{error}</div>}

          <button className="submit-button" type="submit" disabled={!canSubmit}>
            {isSubmitting ? (
              <Loader2 className="spin" size={18} aria-hidden="true" />
            ) : (
              <FlaskConical size={18} aria-hidden="true" />
            )}
            Estimate hinges
            <ArrowRight size={18} aria-hidden="true" />
          </button>
        </form>

        <section className="panel results" aria-live="polite">
          <div className="panel-header">
            <div>
              <span className="section-kicker">Output</span>
              <h2>Results</h2>
            </div>
            {result && <CheckCircle2 className="success-icon" size={22} aria-hidden="true" />}
          </div>

          {!result && !isSubmitting && (
            <div className="empty-state">
              <Activity size={42} aria-hidden="true" />
              <p>Results appear here after estimation.</p>
            </div>
          )}

          {isSubmitting && (
            <div className="empty-state">
              <Loader2 className="spin" size={42} aria-hidden="true" />
              <p>Running estimator...</p>
            </div>
          )}

          {result && (
            <div className="result-list">
              {Object.entries(result.results).map(([key, value]) => (
                <div className="result-row" key={key}>
                  <span>{key}</span>
                  <strong>{value}</strong>
                </div>
              ))}
              {Object.keys(result.results).length === 0 && (
                <pre className="raw-output">{result.rawOutput || 'No structured output returned.'}</pre>
              )}
            </div>
          )}
        </section>
      </section>
    </main>
  );
}

type FileInputProps = {
  id: string;
  label: string;
  file: File | null;
  onChange: (event: ChangeEvent<HTMLInputElement>) => void;
};

function FileInput({ id, label, file, onChange }: FileInputProps) {
  return (
    <label className="file-drop" htmlFor={id}>
      <FileUp size={24} aria-hidden="true" />
      <span>{label}</span>
      <strong>{file?.name ?? 'Choose PDB file'}</strong>
      <input id={id} type="file" accept=".pdb,.ent,.txt" onChange={onChange} required />
    </label>
  );
}

type SegmentedControlProps<T extends string> = {
  value: T;
  options: Array<{ label: string; value: T }>;
  onChange: (value: T) => void;
};

function SegmentedControl<T extends string>({ value, options, onChange }: SegmentedControlProps<T>) {
  return (
    <div className="segmented">
      {options.map((option) => (
        <button
          className={option.value === value ? 'active' : ''}
          key={option.value}
          type="button"
          onClick={() => onChange(option.value)}
        >
          {option.label}
        </button>
      ))}
    </div>
  );
}

export default App;
