# How the Observatory Works — the full loop

From satellite photons to future predictions. **Solid boxes are built & validated;
dashed are in progress; dotted are the roadmap.** The two loops are the point:
validation never stops improving accuracy, and every new season tests the
prediction model against reality.

```mermaid
flowchart TB
    subgraph SENSE["1 · SENSE — the eyes in orbit"]
        S2["Sentinel-2 optical<br/>(area · every few days · cloud-blocked)"]
        S1["Sentinel-1 radar<br/>(sees through cloud & polar night)"]
        IS2["ICESat-2 laser<br/>(pond depth samples)"]
    end

    subgraph DETECT["2 · DETECT — find the water"]
        NDWI["NDWI + shadow test + hysteresis<br/>(published Moussavi method)"]
        CLOUD["Cloud-robust scene selection<br/>(actual-cloud-over-shelf gates)"]
    end

    subgraph VALIDATE["3 · VALIDATE — never trust yourself"]
        BLIND["Blind human labels<br/>(GVI p=0.63 · LarsenC p=0.60)"]
        XSAT["Cross-satellite check<br/>(Landsat-8, r=0.92)"]
        LIT["Literature match<br/>(Banwell 2021 record year)"]
        CORR["Corrections & uncertainty bands<br/>(e.g. slush factor ~×0.6 on Larsen C)"]
    end

    subgraph RECORD["4 · RECORD — the measurement baseline"]
        MULTI["10-shelf, 9-season melt record<br/>(one fixed grid per shelf, comparable years)"]
        PROD["Open data products<br/>(COG + STAC, DOI'd dataset)"]
    end

    subgraph EXPLAIN["5 · EXPLAIN — what drives the melt?"]
        ERA5["Climate records (ERA5)<br/>temperature · wind · sunshine per shelf-season"]
        PATT["Patterns<br/>regional sync (GVI+Wilkins '19-20) ·<br/>new-record shelves ('25-26) · trends"]
    end

    subgraph PREDICT["6 · PREDICT — the payoff"]
        MODEL["Train: melt = f(climate)<br/>on the whole multi-shelf record"]
        FCST["Forecast next season's melt<br/>from climate forecasts"]
        TEST["Every new season TESTS the model<br/>prediction vs observation"]
    end

    subgraph WARN["7 · SHARE & WARN"]
        VULN["Vulnerability signals<br/>% shelf ponded · sudden lake drainage<br/>(hydrofracture precursors)"]
        DASH["Open dashboard + API<br/>for researchers, modelers, anyone"]
    end

    S2 --> NDWI
    S1 -.->|"M3: fusion = gap-free, year-round"| NDWI
    IS2 -.->|"M4: depth → volume km³"| RECORD
    NDWI --> CLOUD --> RECORD
    RECORD --> BLIND & XSAT & LIT
    BLIND & XSAT & LIT --> CORR
    CORR -->|"accuracy loop — runs forever"| DETECT
    CORR --> RECORD
    RECORD --> PROD
    RECORD --> PATT
    ERA5 --> PATT
    PATT --> MODEL --> FCST --> TEST
    TEST -->|"model wrong? retrain"| MODEL
    RECORD -->|"each new season"| TEST
    RECORD --> VULN
    FCST --> VULN
    PROD --> DASH
    VULN --> DASH

    classDef done fill:#1f93cf,stroke:#0e5b86,color:#04121d,font-weight:bold
    classDef doing fill:#a7d9f2,stroke:#1f93cf,color:#0e2a3d,stroke-dasharray: 6 3
    classDef future fill:none,stroke:#8fa3b8,color:#5d7185,stroke-dasharray: 2 3

    class S2,NDWI,CLOUD,BLIND,XSAT,LIT,CORR done
    class MULTI,PROD,RECORD doing
    class S1,IS2,ERA5,PATT,MODEL,FCST,TEST,VULN,DASH future
```

## Where your understanding was spot-on — and the three missing pieces

**You said:** *measure melt → find patterns across shelves → push accuracy up → train a model to predict.* Correct spine.

**Missing piece 1 — weather comes in late.** The satellites don't tell us weather;
they show the water. The climate data (ERA5: temperature, wind, sunshine) joins at
stage 5, where we ask *"what conditions produced the melt we measured?"* That
pairing is the training data for prediction.

**Missing piece 2 — accuracy is a loop, not a gate.** Validation (stage 3) never
finishes. Today's example: blind-labeling Larsen C revealed its "water" is ~40%
slush, so east-side shelves get a correction band. Every such finding flows back
into the detector and the records. The system gets more trustworthy every cycle.

**Missing piece 3 — predictions get graded.** A model that forecasts next summer's
melt is only science if next summer *tests* it (stage 6's feedback arrow). The
observatory keeps observing, so every season is a fresh exam — bad models get
retrained, good ones earn trust.

**Current position:** stages 1–4 built and validated for optical; stage 4's data
products started today (M1). Stage 5 (ERA5) is the next big science step; 6 and 7
stand on everything before them.
