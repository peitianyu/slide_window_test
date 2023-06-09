#ifndef __PATTERN_ORDERING_H__
#define __PATTERN_ORDERING_H__

// for each variable, we need to know its dimension and its index in the H matrix.
struct PatternOrdering
{
public:
    inline PatternOrdering &SetDim(int val)
    {
        m_dim = val;
        return *this;
    }

    inline PatternOrdering &SetIdx(int val)
    {
        m_idx = val;
        return *this;
    }

    inline int Dim() const { return m_dim; }
    inline int Idx() const { return m_idx; }
private:
    int m_dim;
    int m_idx; // The index of this variable in the H matrix.
};



#endif // !__PATTERN_ORDERING_H__