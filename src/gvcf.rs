use crate::{CallError, CalledSite, GvcfSite, Result};

pub struct GvcfBlocker {
    thresholds: Box<[u32]>,
    block: Option<Block>,
    previous: Option<(usize, u64)>,
}

struct Block {
    site: CalledSite,
    end: u64,
    minimum_depth: u32,
    range: usize,
}

impl GvcfBlocker {
    pub fn new(thresholds: impl Into<Box<[u32]>>) -> Result<Self> {
        let thresholds = thresholds.into();
        if thresholds.is_empty() || thresholds.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(CallError::InvalidGvcfThresholds);
        }
        Ok(Self {
            thresholds,
            block: None,
            previous: None,
        })
    }

    pub fn thresholds(&self) -> &[u32] {
        &self.thresholds
    }

    pub fn push(
        &mut self,
        mut site: CalledSite,
        mut emit: impl FnMut(CalledSite) -> Result<()>,
    ) -> Result<()> {
        let coordinate = (site.reference_sequence_id(), site.position());
        if self.previous.is_some_and(|previous| coordinate < previous) {
            return Err(CallError::InvalidGvcfOrder);
        }
        self.previous = Some(coordinate);

        let minimum_depth = site
            .samples()
            .iter()
            .map(|sample| sample.evidence().depth())
            .min()
            .ok_or(CallError::InvalidSampleCount)?;
        let range = self
            .thresholds
            .partition_point(|threshold| *threshold <= minimum_depth);
        let can_collapse = is_reference(&site) && range != 0;
        if can_collapse
            && site.samples().iter().any(|sample| {
                sample
                    .phred_likelihoods()
                    .is_some_and(|values| values.len() != 3)
            })
        {
            return Err(CallError::InvalidGvcfLikelihoods);
        }
        let needs_flush = self.block.as_ref().is_some_and(|block| {
            !can_collapse
                || block.range != range
                || block.site.reference_sequence_id() != site.reference_sequence_id()
                || site.position() > block.end.saturating_add(1)
        });

        if needs_flush {
            let same_end = self.block.as_ref().is_some_and(|block| {
                block.site.reference_sequence_id() == site.reference_sequence_id()
                    && block.end == site.position()
                    && block.end > block.site.position()
            });
            if let Some(block) = &mut self.block
                && same_end
            {
                block.end -= 1;
            }
            if let Some(block) = self.take_block() {
                emit(block)?;
            }
        }

        if can_collapse {
            match &mut self.block {
                Some(block) => block.extend(&site, minimum_depth)?,
                None => {
                    self.block = Some(Block {
                        end: site.position(),
                        site,
                        minimum_depth,
                        range,
                    });
                }
            }
        } else {
            if is_reference(&site) && minimum_depth != 0 {
                site.gvcf = Some(GvcfSite::new(None, minimum_depth, false));
            }
            emit(site)?;
        }
        Ok(())
    }

    pub fn finish(mut self, mut emit: impl FnMut(CalledSite) -> Result<()>) -> Result<()> {
        if let Some(site) = self.take_block() {
            emit(site)?;
        }
        Ok(())
    }

    fn take_block(&mut self) -> Option<CalledSite> {
        self.block.take().map(Block::finish)
    }
}

impl Block {
    fn extend(&mut self, site: &CalledSite, minimum_depth: u32) -> Result<()> {
        if self.site.samples.len() != site.samples.len() {
            return Err(CallError::GvcfSampleCountMismatch);
        }
        self.end = site.position();
        self.minimum_depth = self.minimum_depth.min(minimum_depth);
        for (sample, next) in self.site.samples.iter_mut().zip(site.samples()) {
            let depth = sample.evidence.depth().min(next.evidence().depth());
            sample.evidence.set_depth(depth);
            match (&mut sample.phred_likelihoods, next.phred_likelihoods()) {
                (Some(current), Some(next)) if current.len() == 3 && next.len() == 3 => {
                    if current[1] > next[1] {
                        current[1] = next[1];
                        current[2] = next[2];
                    } else if current[1] == next[1] {
                        current[2] = current[2].min(next[2]);
                    }
                }
                (Some(_), Some(_)) => return Err(CallError::InvalidGvcfLikelihoods),
                (current, None) => *current = None,
                (None, Some(_)) => {}
            }
        }
        Ok(())
    }

    fn finish(mut self) -> CalledSite {
        let end_position = (self.end > self.site.position()).then_some(self.end);
        self.site.quality = None;
        self.site.indel_summary = None;
        self.site.annotations = None;
        self.site.gvcf = Some(GvcfSite::new(end_position, self.minimum_depth, true));
        self.site
    }
}

fn is_reference(site: &CalledSite) -> bool {
    !site.is_variant()
        && site.alternates().is_empty()
        && site.samples().iter().all(|sample| {
            sample
                .genotype()
                .is_none_or(|genotype| genotype.iter().all(|allele| *allele == 0))
        })
}

#[cfg(test)]
mod tests;
