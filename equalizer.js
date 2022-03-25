/*
Biquad algorithms are taken from:
https://github.com/jaakkopasanen/AutoEq/blob/master/biquad.py
https://github.com/mohayonao/biquad-coeffs/tree/master/packages/biquad-coeffs-cookbook
*/

Equalizer = (function() {
    let config = {
        DefaultSampleRate: 48000
    };

    let lowshelf = function (freq, q, gain, sampleRate) {
        freq = freq / (sampleRate || config.DefaultSampleRate);
        freq = Math.max(1e-6, Math.min(freq, 1));
        q    = Math.max(1e-4, Math.min(q, 1000));
        gain = Math.max(-40, Math.min(gain, 40));

        let w0 = 2 * Math.PI * freq;
        let sin = Math.sin(w0);
        let cos = Math.cos(w0);
        let a = Math.pow(10, (gain / 40));
        let alpha = sin / (2 * q);
        let alphamod = (2 * Math.sqrt(a) * alpha) || 0;

        let a0 =          ((a+1) + (a-1) * cos + alphamod);
        let a1 = -2 *     ((a-1) + (a+1) * cos           );
        let a2 =          ((a+1) + (a-1) * cos - alphamod);
        let b0 =      a * ((a+1) - (a-1) * cos + alphamod);
        let b1 =  2 * a * ((a-1) - (a+1) * cos           );
        let b2 =      a * ((a+1) - (a-1) * cos - alphamod);

        return [ 1.0, a1/a0, a2/a0, b0/a0, b1/a0, b2/a0 ];
    };

    let highshelf = function (freq, q, gain, sampleRate) {
        freq = freq / (sampleRate || config.DefaultSampleRate);
        freq = Math.max(1e-6, Math.min(freq, 1));
        q    = Math.max(1e-4, Math.min(q, 1000));
        gain = Math.max(-40, Math.min(gain, 40));

        let w0 = 2 * Math.PI * freq;
        let sin = Math.sin(w0);
        let cos = Math.cos(w0);
        let a = Math.pow(10, (gain / 40));
        let alpha = sin / (2 * q);
        let alphamod = (2 * Math.sqrt(a) * alpha) || 0;

        let a0 =          ((a+1) - (a-1) * cos + alphamod);
        let a1 =  2 *     ((a-1) - (a+1) * cos           );
        let a2 =          ((a+1) - (a-1) * cos - alphamod);
        let b0 =      a * ((a+1) + (a-1) * cos + alphamod);
        let b1 = -2 * a * ((a-1) + (a+1) * cos           );
        let b2 =      a * ((a+1) + (a-1) * cos - alphamod);

        return [ 1.0, a1/a0, a2/a0, b0/a0, b1/a0, b2/a0 ];
    };

    let peaking = function (freq, q, gain, sampleRate) {
        freq = freq / (sampleRate || config.DefaultSampleRate);
        freq = Math.max(1e-6, Math.min(freq, 1));
        q    = Math.max(1e-4, Math.min(q, 1000));
        gain = Math.max(-40, Math.min(gain, 40));

        let w0 = 2 * Math.PI * freq;
        let sin = Math.sin(w0);
        let cos = Math.cos(w0);
        let a = Math.pow(10, (gain / 40));
        let alpha = sin / (2 * q);

        let a0 =  1 + alpha / a;
        let a1 = -2 * cos;
        let a2 =  1 - alpha / a;
        let b0 =  1 + alpha * a;
        let b1 = -2 * cos;
        let b2 =  1 - alpha * a;

        return [ 1.0, a1/a0, a2/a0, b0/a0, b1/a0, b2/a0 ];
    };

    let calc_gains = function (freqs, coeffs, sampleRate) {
        sampleRate = sampleRate || config.DefaultSampleRate;
        let gains = new Array(freqs.length).fill(0);

        for (let i = 0; i < coeffs.length; ++i) {
            let [ a0, a1, a2, b0, b1, b2] = coeffs[i];
            for (let j = 0; j < freqs.length; ++j) {
                let w = 2 * Math.PI * freqs[j] / sampleRate;
                let phi = 4 * Math.pow(Math.sin(w / 2), 2);
                let c = (
                    10 * Math.log10(Math.pow(b0 + b1 + b2, 2) +
                        (b0 * b2 * phi - (b1 * (b0 + b2) + 4 * b0 * b2)) * phi) -
                    10 * Math.log10(Math.pow(a0 + a1 + a2, 2) +
                        (a0 * a2 * phi - (a1 * (a0 + a2) + 4 * a0 * a2)) * phi));
                gains[j] += c;
            }
        }
        return gains;
    };

    let calc_preamp = function (fr1, fr2) {
        let max1 = -Infinity;
        let max2 = -Infinity;
        for (let i = 0; i < fr1.length; ++i) {
            max1 = Math.max(max1, fr1[i][1]);
        }
        for (let i = 0; i < fr2.length; ++i) {
            max2 = Math.max(max2, fr2[i][1]);
        }
        return max1 - max2;
    };

    let apply = function (fr, filters, sampleRate) {
        let freqs = new Array(fr.length).fill(null);
        for (let i = 0; i < fr.length; ++i) {
            freqs[i] = fr[i][0];
        }
        let coeffs = filters.map(f => {
            if (!f.freq || !f.gain || !f.q) {
                return null;
            } else if (f.type === "LSQ") {
                return lowshelf(f.freq, f.q, f.gain, sampleRate);
            } else if (f.type === "HSQ") {
                return highshelf(f.freq, f.q, f.gain, sampleRate);
            } else if (f.type === "PK") {
                return peaking(f.freq, f.q, f.gain, sampleRate);
            }
            return null;
        }).filter(f => f);
        let gains = calc_gains(freqs, coeffs, sampleRate);
        let fr_eq = new Array(fr.length).fill(null);
        for (let i = 0; i < fr.length; ++i) {
            fr_eq[i] = [fr[i][0], fr[i][1] + gains[i]];
        }
        return fr_eq;
    };

    return {
        config,
        lowshelf,
        highshelf,
        peaking,
        calc_gains,
        calc_preamp,
        apply
    }
})();

