# GRL reviewer analyses

This folder contains analyses added in response to reviewer comments on the GRL manuscript. It is kept separate from the submitted-manuscript workflow so that each new analysis can be developed and documented without changing the archived analysis code.

## Structure

- `00_Data/`: descriptions of required raw and processed inputs. Large datasets are not stored in Git.
- `01_Codes/`: R scripts used for reviewer analyses.
- `02_Results/`: figures and tables generated for the response and revised manuscript.
- `03_Documentation/`: reviewer-comment tracking and data requirements.

## Working procedure

For each reviewer comment:

1. Record the comment and proposed response in `03_Documentation/REVIEWER_COMMENTS.md`.
2. Identify the necessary input data before writing or changing analysis code.
3. Add the analysis to `01_Codes/` while preserving the existing code style.
4. Save compact response-ready outputs under `02_Results/`.
5. Record the result and its interpretation beside the corresponding reviewer comment.

The current scripts retain their original OSC paths and code annotations. Paths can be updated when the working data location is confirmed.
