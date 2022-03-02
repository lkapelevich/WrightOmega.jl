# approximation algorithm from ``Algorithm 917: Complex Double-Precision
# Evaluation of the Wright ω Function'' by Piers W. Laurence, Robert M. Corless
# and David J. Jeffrey
# Series approximations can be found in the book chapter ``The Wright ω Function''
# by Corless, R. M. and Jeffrey, D. J.
#
# see also https://github.com/emsr/tr29124_test

module WrightOmega

export wrightomega

# TODO handle infs
# TODO don't use two methods
# FIXME test failures around zero and regularization

function wrightomega(z::T) where {T <: Real}
    if z < -2
        t = exp(z)
        w = t * (1 + t * (-1 + t * (T(3) / 2 + t * (-T(8) / 3 + T(125) / 24 * t))))
    elseif z < 1 + T(π)
        z1 = z - 1
        w = 1 + z1 / 2 * (1 + z1 / 8 * (1 + z1 / 12 * (-1 + z1 / 16 *
            (-1 + z1 * T(13) / 20))))
    else
        l = log(z)
        w = z + l * (-1 + (1 + (l / 2 - 1 + (l * (l / 3 - T(3) / 2) + 1) / z) / z) / z)
    end
    # refinement steps
    for _ in 1:5
        r = z - w - log(w)
        w1 = w + 1
        t = w1 * (w1 + T(2) / 3 * r)
        w *= 1 + r / w1 * (t - r / 2) / (t - r)
        fscn = abs(r^4 * (2 * w * (w - 4) - 1))
        if fscn < eps(float(T)) * 72 * w1^6
            break
        end
    end
    return w::float(T)
end

function wrightomega(z::Complex{T}) where {T <: Real}
    x = real(z)
    y = imag(z)
    Tπ = float(T)(π)
    if x == -1 && abs(y) == Tπ
        return zero(Complex{T})
    end
    # region 1
    if -2 < x <= 1 && 1 < y < 2Tπ
        imp = im * conj(sqrt(2 * conj(z + 1 - im * Tπ)))
        w = -1 + imp * (1 + imp / 3 * (-1 + imp / 6 * (inv(T(2)) + imp / 15 *
            (1 + imp / 16))))
    # region 2
    elseif -2 < x <= 1 && -2Tπ < y < -1
        imp = im * conj(sqrt(2 * conj(z + 1 - im * Tπ)))
        w = -1 + imp * (-1 + imp / 3 * (-1 + imp / 6 * (-inv(T(2)) + imp / 15 *
            (1 - imp / 16))))
    # region 3
    elseif x < -2 && -Tπ < y <= Tπ
        t = exp(z)
        w = t * (1 + t * (-1 + t * (T(3) / 2 + t * (-T(8) / 3 + T(125) / 24 * t))))
    # region 4
    elseif (-2 < x <= 1 && -1 <= y <= 1) || (-1 < x && (x - 1)^2 + y^2 <= Tπ^2)
        z1 = z - 1
        w = 1 + z1 / 2 * (1 + z1 / 8 * (1 + z1 / 12 * (-1 + z1 / 16 *
            (-1 + z1 * T(13) / 20))))
    # region 5
    elseif x <= -2 && Tπ < y && y - Tπ <= -3 / T(4) * (x + 1)
        t = z - im * Tπ
        l = log(-t)
        w = t - l + l * (1 + (l / 2 - 1 + (l^2 / 3 - 3 * l / 2 + 1) / t) / t) / t
    # region 6
    elseif x <= -2 && y <= -Tπ && y + Tπ >= T(3) / 4 * (x + 1)
        t = z + im * Tπ
        l = log(-t)
        w = t - l + l * (1 + (l / 2 - 1 + (l^2 / 3 - 3 * l / 2 + 1) / t) / t) / t
    # region 7
    else
        l = log(z)
        w = z + l * (-1 + (1 + (l / 2 - 1 + (l * (l / 3 - T(3) / 2) + 1) / z) / z) / z)
    end
    if x <= -0.99 && (abs(y - Tπ) <= 0.01 || abs(y + Tπ) <= 0.01)
        s = -1
        if abs(y - Tπ) <= 0.01
            # TODO setrounding doesn't work for non-BigFloat, decide what to do
            # T == BigFloat && setrounding(T, RoundUp)
            z = x + im * (y - Tπ)
        else
            z = x + im * (y + Tπ)
        end
    else
        s = 1
    end
    w *= s
    # refinement steps
    for _ in 1:5
        r = z - s * w - log(Complex(w))
        w1 = s * w + 1
        t = w1 * (w1 + T(2) / 3 * r)
        w *= 1 + r / w1 * (t - r / 2) / (t - r)
        fscn = abs(r^4 * (2 * w * (w - 4) - 1))
        if fscn < eps(float(T)) * 72 * abs(w1)^6
            break
        end
    end
    w *= s
    return w::Complex{float(T)}
end

end
