---
title: Help
description: 'Where to ask a question or report a bug, how to contribute, what is planned, and the code of conduct.'
---

## Where to ask

Ask on [GitHub issues](https://github.com/slamballais/QDECR/issues), for questions as well as bugs. Search the open and closed issues first: most problems have come up before, and so have their fixes.

Many problems that look like QDECR errors come from FreeSurfer's setup or from the data frame. [Troubleshooting and FAQ](/tutorials/troubleshooting) collects the common ones, with the error messages they produce.

QDECR is maintained in the authors' own time. Issues are answered as soon as we can, but not on a schedule.

## Reporting a bug

A report we can act on says what you ran, what happened, and on what. Please include:

- The call you ran, such as `qdecr_fastlm(...)` with all its arguments.
- The full error message or the unexpected output, as text rather than a screenshot.
- The result printed with `print(out)`, if the analysis got that far.
- Your R session and the FreeSurfer build:

```r
sessionInfo()
packageVersion("QDECR")
```

```bash
cat "$FREESURFER_HOME/build-stamp.txt"
```

- Your operating system, and whether R runs inside WSL2, a container or on a cluster.
- The number of subjects and imputed datasets, and the `n_cores` you used.

Leave out anything that identifies participants. If the problem needs data to reproduce, a handful of subjects from a public dataset is ideal.

## Contributing

Contributions are welcome, from a fixed typo to a new model.

- **Documentation.** Every page on this site has an "Edit this page" link in its footer that opens its source on GitHub. The site lives in `website/` on the `master` branch, next to the release it documents, so site changes branch off `master`. The words it uses for things are those of the [glossary](/glossary).
- **Bugs and features.** Open an issue first for anything larger than a small fix, so we can agree on the approach before you write the code. Package changes branch off `dev`, where development happens, and are merged into `master` at a release.
- **Pull requests.** Keep each one small and about one thing: a pull request that does one thing is reviewed quickly, one that does five things waits. Expect a round or two of review.

## Roadmap

QDECR is stable: version 0.9.0 is the current release, and 0.10.0 is in development. The status of what has been asked for over the years:

| Feature | Status |
|---|---|
| Weighted regression | Done in 0.9.0: the `weights` argument of `qdecr_fastlm`. |
| Your own surface maps | Done in 0.8.5: the `custom_measure` argument. |
| Package upkeep | In development for 0.10.0: citation file, licence file, automated checks and tests. |
| Targets other than `fsaverage` | A bug fix planned for 0.10.0. The `target` argument is not yet used everywhere it should be. |
| Linear mixed models | Not planned in QDECR. [verywise](https://github.com/SereDef/verywise) fits vertex-wise linear mixed models for longitudinal and multi-site data, with the same FreeSurfer correction. |
| Permutation testing | Not planned. |
| Data from other software (CIVET, ANTs) | Not planned. Any map that is resampled to `fsaverage` and saved as an MGH file per subject can already be analysed with `custom_measure`. |
| Windows without WSL2 | Not planned: QDECR runs its parallel work in forked processes, which Windows does not have. [Use WSL2](/get-started#windows-through-wsl2). |

## Code of conduct

Everyone who takes part in QDECR's issues, pull requests and other spaces agrees to the code below: the [Contributor Covenant](https://www.contributor-covenant.org/), version 3.0. The project's maintainers are its Community Moderators.

### Our pledge

We pledge to make our community welcoming, safe, and equitable for all.

We are committed to fostering an environment that respects and promotes the dignity, rights, and contributions of all individuals, regardless of characteristics including race, ethnicity, caste, color, age, physical characteristics, neurodiversity, disability, sex or gender, gender identity or expression, sexual orientation, language, philosophy or religion, national or social origin, socio-economic position, level of education, or other status. The same privileges of participation are extended to everyone who participates in good faith and in accordance with this Covenant.

### Encouraged behaviors

While acknowledging differences in social norms, we all strive to meet our community's expectations for positive behavior. We also understand that our words and actions may be interpreted differently than we intend based on culture, background, or native language.

With these considerations in mind, we agree to behave mindfully toward each other and act in ways that center our shared values, including:

1. Respecting the **purpose of our community**, our activities, and our ways of gathering.
2. Engaging **kindly and honestly** with others.
3. Respecting **different viewpoints** and experiences.
4. **Taking responsibility** for our actions and contributions.
5. Gracefully giving and accepting **constructive feedback**.
6. Committing to **repairing harm** when it occurs.
7. Behaving in other ways that promote and sustain the **well-being of our community**.

### Restricted behaviors

We agree to restrict the following behaviors in our community. Instances, threats, and promotion of these behaviors are violations of this Code of Conduct.

1. **Harassment.** Violating explicitly expressed boundaries or engaging in unnecessary personal attention after any clear request to stop.
2. **Character attacks.** Making insulting, demeaning, or pejorative comments directed at a community member or group of people.
3. **Stereotyping or discrimination.** Characterizing anyone's personality or behavior on the basis of immutable identities or traits.
4. **Sexualization.** Behaving in a way that would generally be considered inappropriately intimate in the context or purpose of the community.
5. **Violating confidentiality**. Sharing or acting on someone's personal or private information without their permission.
6. **Endangerment.** Causing, encouraging, or threatening violence or other harm toward any person or group.
7. Behaving in other ways that **threaten the well-being** of our community.

#### Other restrictions

1. **Misleading identity.** Impersonating someone else for any reason, or pretending to be someone else to evade enforcement actions.
2. **Failing to credit sources.** Not properly crediting the sources of content you contribute.
3. **Promotional materials**. Sharing marketing or other commercial content in a way that is outside the norms of the community.
4. **Irresponsible communication.** Failing to responsibly present content which includes, links or describes any other restricted behaviors.

### Reporting an issue

Tensions can occur between community members even when they are trying their best to collaborate. Not every conflict represents a code of conduct violation, and this Code of Conduct reinforces encouraged behaviors and norms that can help avoid conflicts and minimize harm.

When an incident does occur, it is important to report it promptly. To report a possible violation, **email [contact@qdecr.com](mailto:contact@qdecr.com). Only the maintainers read it.**

Community Moderators take reports of violations seriously and will make every effort to respond in a timely manner. They will investigate all reports of code of conduct violations, reviewing messages, logs, and recordings, or interviewing witnesses and other participants. Community Moderators will keep investigation and enforcement actions as transparent as possible while prioritizing safety and confidentiality. In order to honor these values, enforcement actions are carried out in private with the involved parties, but communicating to the whole community may be part of a mutually agreed upon resolution.

### Addressing and repairing harm

If an investigation by the Community Moderators finds that this Code of Conduct has been violated, the following enforcement ladder may be used to determine how best to repair harm, based on the incident's impact on the individuals involved and the community as a whole. Depending on the severity of a violation, lower rungs on the ladder may be skipped.

1. Warning
   1. Event: A violation involving a single incident or series of incidents.
   2. Consequence: A private, written warning from the Community Moderators.
   3. Repair: Examples of repair include a private written apology, acknowledgement of responsibility, and seeking clarification on expectations.
2. Temporarily Limited Activities
   1. Event: A repeated incidence of a violation that previously resulted in a warning, or the first incidence of a more serious violation.
   2. Consequence: A private, written warning with a time-limited cooldown period designed to underscore the seriousness of the situation and give the community members involved time to process the incident. The cooldown period may be limited to particular communication channels or interactions with particular community members.
   3. Repair: Examples of repair may include making an apology, using the cooldown period to reflect on actions and impact, and being thoughtful about re-entering community spaces after the period is over.
3. Temporary Suspension
   1. Event: A pattern of repeated violation which the Community Moderators have tried to address with warnings, or a single serious violation.
   2. Consequence: A private written warning with conditions for return from suspension. In general, temporary suspensions give the person being suspended time to reflect upon their behavior and possible corrective actions.
   3. Repair: Examples of repair include respecting the spirit of the suspension, meeting the specified conditions for return, and being thoughtful about how to reintegrate with the community when the suspension is lifted.
4. Permanent Ban
   1. Event: A pattern of repeated code of conduct violations that other steps on the ladder have failed to resolve, or a violation so serious that the Community Moderators determine there is no way to keep the community safe with this person as a member.
   2. Consequence: Access to all community spaces, tools, and communication channels is removed. In general, permanent bans should be rarely used, should have strong reasoning behind them, and should only be resorted to if working through other remedies has failed to change the behavior.
   3. Repair: There is no possible repair in cases of this severity.

This enforcement ladder is intended as a guideline. It does not limit the ability of Community Managers to use their discretion and judgment, in keeping with the best interests of our community.

### Scope

This Code of Conduct applies within all community spaces, and also applies when an individual is officially representing the community in public or other spaces. Examples of representing our community include using an official email address, posting via an official social media account, or acting as an appointed representative at an online or offline event.

### Attribution

This Code of Conduct is adapted from the Contributor Covenant, version 3.0, permanently available at [https://www.contributor-covenant.org/version/3/0/](https://www.contributor-covenant.org/version/3/0/).

Contributor Covenant is stewarded by the Organization for Ethical Source and licensed under CC BY-SA 4.0. To view a copy of this license, visit [https://creativecommons.org/licenses/by-sa/4.0/](https://creativecommons.org/licenses/by-sa/4.0/)

For answers to common questions about Contributor Covenant, see the FAQ at [https://www.contributor-covenant.org/faq](https://www.contributor-covenant.org/faq). Translations are provided at [https://www.contributor-covenant.org/translations](https://www.contributor-covenant.org/translations). Additional enforcement and community guideline resources can be found at [https://www.contributor-covenant.org/resources](https://www.contributor-covenant.org/resources). The enforcement ladder was inspired by the work of [Mozilla's code of conduct team](https://github.com/mozilla/inclusion).
